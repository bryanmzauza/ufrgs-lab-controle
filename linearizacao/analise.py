"""Linearização da bancada de seis tanques e sintonia das quatro malhas.

1. O OpenModelica lineariza a SixTanks_MA (a SixTanks em malha aberta, parada no
   ponto de operação do MF_ST) e devolve A, B, C, D.
2. Funções de transferência da planta, RGA e zeros de transmissão.
3. Sintonia SIMC (PI): C3/C4 primeiro (T3/T4 são integradores), depois C1/C2
   com C3/C4 fechados.
4. Validação: o MF_ST não linear, já em regime, recebe os degraus de SP e é
   comparado com o modelo linear em malha fechada.

Uso:
    python analise.py              # relatório no terminal e validacao.png nesta pasta
    python analise.py --sem-figura

Requer numpy e scipy (matplotlib só para a figura) e o OpenModelica.
"""

import argparse
import csv
import importlib.util
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy import linalg, optimize, signal

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "simulador"))
from servidor import achar_omc  # noqa: E402

ARQUIVO_MO = os.path.normpath(os.path.join(AQUI, "..", "modelos", "modelo_6tanques_v12.mo"))
PASTA = os.path.join(tempfile.gettempdir(), "seis_tanques_lin")

TAU_C_T3T4 = 10.0   # s, constante de tempo em malha fechada pedida para C3/C4 (integradores sem atraso)
T_DEGRAU = 3000.0   # s, instante dos degraus de SP na validação (o MF_ST já está em regime)
T_VALIDACAO = 1500.0  # s simulados depois do degrau

NIVEIS = ["h1", "h2", "h3", "h4", "h5", "h6"]
VALVULAS = ["FV-01", "FV-02", "FV-03", "FV-04"]


# ---------------------------------------------------------------- OpenModelica
def omc(script, nome):
    exe = achar_omc()
    if not exe:
        sys.exit("OpenModelica (omc) não encontrado.")
    os.makedirs(PASTA, exist_ok=True)
    with open(os.path.join(PASTA, nome + ".mos"), "w", encoding="utf-8") as f:
        f.write(script)
    r = subprocess.run([exe, nome + ".mos"], cwd=PASTA, capture_output=True, text=True, timeout=900)
    if r.returncode != 0 or "Error" in r.stdout:
        sys.exit(f"Falha no OpenModelica ({nome}):\n{r.stdout}\n{r.stderr}")
    return r.stdout


def ler_csv(caminho):
    with open(caminho, newline="") as f:
        linhas = list(csv.reader(f))
    cab, dados = linhas[0], np.array(linhas[1:], dtype=float)
    return {nome: dados[:, i] for i, nome in enumerate(cab)}


def linearizar():
    mo = ARQUIVO_MO.replace("\\", "/")
    omc(f'loadFile("{mo}"); getErrorString();\n'
        'setCommandLineOptions("--linearizationDumpLanguage=python"); getErrorString();\n'
        'linearize(SixTanks_MA, stopTime=1, outputFormat="csv"); getErrorString();\n', "linearizar")
    spec = importlib.util.spec_from_file_location("linearized_model", os.path.join(PASTA, "linearized_model.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    n, m, p, x0, u0, A, B, C, D, sv, iv, ov = mod.linearized_model()
    assert sv == ["P_h1", "P_h2", "P_h3", "P_h4", "P_h5", "P_h6"] and iv == ["u1", "u2", "u3", "u4"], (sv, iv)
    limpa = lambda M: np.where(np.abs(np.array(M, dtype=float)) < 1e-7, 0.0, np.array(M, dtype=float))
    A, B, C = limpa(A), limpa(B), limpa(C)
    res = ler_csv(os.path.join(PASTA, "SixTanks_MA_res.csv"))
    x_op = np.array([res[f"P.x{i}"][-1] for i in range(1, 5)])
    return A, B, C, np.array(x0), x_op


# ---------------------------------------------------------------- ferramentas lineares
def fecha_malhas(A, B, C, malhas, livres):
    """Fecha PIs (forma do PID_ISA com b = 1, Td = 0): u = Kp*(r - y) + Ui, Ti*Ui' = Kp*(r - y).
    malhas: [(saida, entrada, Kp, Ti)]. Entradas do sistema resultante: SPs das malhas e depois as
    entradas livres; saídas: os estados da planta e depois as MVs das malhas."""
    n, L = A.shape[0], len(malhas)
    Sy = np.zeros((L, C.shape[0]))
    Su = np.zeros((B.shape[1], L))
    for l, (yo, ui, _, _) in enumerate(malhas):
        Sy[l, yo] = 1
        Su[ui, l] = 1
    Kc = np.diag([k for _, _, k, _ in malhas])
    Ki = np.diag([k / ti for _, _, k, ti in malhas])
    Cy = Sy @ C
    Acl = np.block([[A - B @ Su @ Kc @ Cy, B @ Su], [-Ki @ Cy, np.zeros((L, L))]])
    Br = np.vstack([B @ Su @ Kc, Ki])
    Bv = np.vstack([B[:, livres], np.zeros((L, len(livres)))])
    Cx = np.hstack([np.eye(n), np.zeros((n, L))])
    Cu = np.hstack([-Kc @ Cy, np.eye(L)])
    Dr = np.vstack([np.zeros((n, L)), Kc])
    return Acl, np.hstack([Br, Bv]), np.vstack([Cx, Cu]), np.hstack([Dr, np.zeros((n + L, len(livres)))])


def ganho_estatico(A, B, C, D):
    return D - C @ np.linalg.solve(A, B)


def zeros_transmissao(A, B, C, D):
    n = A.shape[0]
    M = np.block([[A, B], [C, D]])
    N = np.block([[np.eye(n), np.zeros_like(B)], [np.zeros_like(C), np.zeros_like(D)]])
    z = linalg.eigvals(M, N)
    return np.sort_complex(z[np.isfinite(z) & (np.abs(z) < 1e6)])


def fmt_tf(A, B, C, i, j):
    """G(s) de u_j para y_i em forma de constantes de tempo, sem os polos que o canal não enxerga."""
    z, p, k = signal.ss2zpk(A, B[:, [j]], C[[i], :], np.zeros((1, 1)))
    z, p = list(np.atleast_1d(z[0] if np.ndim(z) > 1 else z)), list(p)
    for zz in list(z):  # cancela polo-zero coincidentes
        d = [abs(zz - pp) for pp in p]
        if d and min(d) < 1e-6 * max(1.0, abs(zz)):
            p.pop(int(np.argmin(d)))
            z.remove(zz)
    if abs(k) < 1e-12:
        return None
    integ = sum(1 for pp in p if abs(pp) < 1e-9)
    p = [pp for pp in p if abs(pp) >= 1e-9]
    K = k * np.prod([-zz for zz in z]).real / np.prod([-pp for pp in p]).real if p else k
    num = "".join(f"({-1 / zz.real:.1f}s + 1)" for zz in z)
    den = ("s" if integ else "") + "".join(f"({-1 / pp.real:.1f}s + 1)" for pp in sorted(p, key=lambda q: q.real))
    return f"{K:+.4g}{' ' + num if num else ''} / {den}"


def resposta_degrau(A, B, C, D, j, t):
    _, y, _ = signal.lsim(signal.StateSpace(A, B[:, [j]], C, D[:, [j]]), np.ones_like(t), t)
    return y


def ajusta_sopdt(t, y):
    """k e^(-θs) / ((τ1 s + 1)(τ2 s + 1)) por mínimos quadrados."""
    def modelo(q):
        k, t1, t2, th = q
        tt = np.maximum(t - th, 0)
        if abs(t1 - t2) < 1e-3:
            t2 = t1 - 1e-3
        return k * (1 - (t1 * np.exp(-tt / t1) - t2 * np.exp(-tt / t2)) / (t1 - t2))
    q0 = [y[-1], 30, 10, 1]
    r = optimize.least_squares(lambda q: modelo(q) - y, q0, bounds=([-np.inf, 0.1, 0.1, 0], [np.inf, 1e4, 1e4, 200]))
    k, t1, t2, th = r.x
    t1, t2 = max(t1, t2), min(t1, t2)
    return k, t1, t2, th, np.max(np.abs(modelo(r.x) - y)) / abs(k)


def simc_pi(k, t1, t2, th):
    """Regra da metade (SOPDT -> FOPDT) e SIMC com τc = θ."""
    tau, theta = t1 + t2 / 2, th + t2 / 2
    tc = theta
    return tau / (k * (tc + theta)), min(tau, 4 * (tc + theta)), tau, theta


def arred(x, sig=2):
    return float(f"{x:.{sig}g}")


# ---------------------------------------------------------------- análise
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sem-figura", action="store_true")
    a = ap.parse_args()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")  # console do Windows (cp1252) não tem τ, θ, λ
    print(f"modelo: {os.path.relpath(ARQUIVO_MO, os.path.dirname(AQUI))} (SixTanks_MA)\n")

    A, B, C, h_op, x_op = linearizar()
    n = A.shape[0]

    print("1. Ponto de operação (equilíbrio da SixTanks_MA)")
    print("   níveis (cm):  " + "  ".join(f"{nm} = {v:.2f}" for nm, v in zip(NIVEIS, h_op)))
    print("   aberturas:    " + "  ".join(f"{nm} = {v:.3f}" for nm, v in zip(VALVULAS, x_op)))
    baixos = [nm for nm, v in zip(NIVEIS, h_op) if v < 5]
    print("   todos os tanques >= 5 cm" if not baixos else f"   ABAIXO DE 5 cm: {', '.join(baixos)}")

    print("\n2. Modelo linear  dh/dt = A dh + B dx   (dh em cm, dx em fração de abertura)")
    np.set_printoptions(precision=4, suppress=True, linewidth=120)
    print("   A =\n" + "\n".join("     " + ln for ln in str(A).splitlines()))
    print("   B =\n" + "\n".join("     " + ln for ln in str(B).splitlines()))
    print("   polos: " + ", ".join(f"{p.real:.4f}" for p in np.sort(np.linalg.eigvals(A).real)))
    print("   T3 e T4 têm polo em zero: a vazão de B-03/B-04 não depende do nível (acima de h_succao_min).")
    print("\n   Funções de transferência G(s) = dh/dx (cm por unidade de abertura, τ em s):")
    for j in range(4):
        for i in range(n):
            tf = fmt_tf(A, B, C, i, j)
            if tf:
                print(f"     {NIVEIS[i]} / {VALVULAS[j]}:  {tf}")

    # --- C3/C4: G = k/s, SIMC para integrador sem atraso: Kc = 1/(k τc), Ti = 4 τc
    print(f"\n3. C3 e C4 (T3/T4 integradores, sem atraso), SIMC com τc = {TAU_C_T3T4:.0f} s")
    sint = {}
    for c, (i, j) in (("C3", (2, 2)), ("C4", (3, 3))):
        k = B[i, j]
        kp, ti = arred(1 / (k * TAU_C_T3T4)), arred(4 * TAU_C_T3T4)
        sint[c] = (i, j, kp, ti)
        print(f"   {c}: {NIVEIS[i]}/{VALVULAS[j]} = {k:.4f}/s  ->  Kp = {kp}, Ti = {ti:.0f} s")
    internas = [sint["C3"], sint["C4"]]

    # --- h1, h2 com C3/C4 fechados: ganho estático, RGA e zeros
    Ai, Bi, Ci, Di = fecha_malhas(A, B, C, internas, livres=[0, 1])
    idx_u, idx_y = [2, 3], [0, 1]  # entradas x1, x2 (depois dos 2 SPs); saídas h1, h2
    K = ganho_estatico(Ai, Bi[:, idx_u], Ci[idx_y], Di[np.ix_(idx_y, idx_u)])
    rga = K * np.linalg.inv(K).T
    z = zeros_transmissao(Ai, Bi[:, idx_u], Ci[idx_y], Di[np.ix_(idx_y, idx_u)])
    print("\n4. h1, h2 em relação a x1, x2 com C3 e C4 fechados")
    print(f"   ganho estático K (cm por unidade de abertura) = [[{K[0,0]:.2f}, {K[0,1]:.2f}], [{K[1,0]:.2f}, {K[1,1]:.2f}]]")
    print(f"   RGA: λ11 = {rga[0,0]:.3f} (h1-FV-01, h2-FV-02),  λ12 = {rga[0,1]:.3f} (h1-FV-02, h2-FV-01)")
    cruzado = rga[0, 1] > rga[0, 0]
    print(f"   pareamento indicado: {'cruzado (C1 -> FV-02, C2 -> FV-01)' if cruzado else 'direto (C1 -> FV-01, C2 -> FV-02)'}")
    zs = ", ".join(f"{v.real:+.4f}" + (f"{v.imag:+.4f}j" if abs(v.imag) > 1e-9 else "") for v in z)
    print(f"   zeros de transmissão (1/s): {zs}")
    spd = [v for v in z if v.real > 1e-9]
    if spd:
        print(f"   zero no semiplano direito em {min(v.real for v in spd):.3f} 1/s: limita a banda passante conjunta de C1/C2 "
              f"a ~{min(v.real for v in spd) / 2:.3f} rad/s (r1 + r2 < 1)")

    # --- C1/C2: SISO com C3/C4 fechados e a outra malha aberta
    print("\n5. C1 e C2 (SISO com C3/C4 fechados e a outra malha aberta): ajuste SOPDT + SIMC PI (τc = θ)")
    t = np.linspace(0, 1500, 3001)
    pares = {"C1": (0, 1 if cruzado else 0), "C2": (1, 0 if cruzado else 1)}
    for c, (i, j) in pares.items():
        y = resposta_degrau(Ai, Bi, Ci, Di, 2 + j, t)[:, i]
        k, t1, t2, th, erro = ajusta_sopdt(t, y)
        kc, ti, tau, theta = simc_pi(k, t1, t2, th)
        sint[c] = (i, j, arred(kc), arred(ti))
        print(f"   {c}: {NIVEIS[i]}/{VALVULAS[j]} ≈ {k:.2f} e^(-{th:.1f}s) / (({t1:.1f}s + 1)({t2:.1f}s + 1))"
              f"  [erro máx. {100 * erro:.1f} %]")
        print(f"       metade: τ = {tau:.1f} s, θ = {theta:.1f} s  ->  Kp = {sint[c][2]}, Ti = {sint[c][3]:.0f} s")

    # --- malha fechada completa
    malhas = [sint[c] for c in ("C1", "C2", "C3", "C4")]
    Af, Bf, Cf, Df = fecha_malhas(A, B, C, malhas, livres=[])
    pf = np.linalg.eigvals(Af)
    print(f"\n6. Malha fechada linear com as quatro malhas: polo mais lento {max(pf.real):.4f} 1/s "
          f"({'estável' if max(pf.real) < 0 else 'INSTÁVEL'}), menor amortecimento "
          f"{min(-p.real / abs(p) for p in pf if abs(p) > 1e-12):.2f}")

    print("\n   Para o MF_ST:")
    for c in ("C1", "C2", "C3", "C4"):
        print(f"     {c}: Kp = {sint[c][2]}, Ti = {sint[c][3]:.0f}, Td = 0")

    validar(A, B, C, h_op, x_op, sint, Af, Bf, Cf, Df, a.sem_figura)


def validar(A, B, C, h_op, x_op, sint, Af, Bf, Cf, Df, sem_figura):
    """MF_ST não linear com a sintonia acima, degraus de SP em T_DEGRAU, contra o modelo linear."""
    mo = ARQUIVO_MO.replace("\\", "/")
    ov = [f"TimeYset={T_DEGRAU}"] + [f"{c}.{p}={v}" for c, (_, _, kp, ti) in sint.items() for p, v in (("Kp", kp), ("Ti", ti))]
    omc(f'loadFile("{mo}"); getErrorString();\n'
        f'simulate(MF_ST, stopTime={T_DEGRAU + T_VALIDACAO}, numberOfIntervals={int((T_DEGRAU + T_VALIDACAO) / 2)}, '
        f'outputFormat="csv", simflags="-override={",".join(ov)}"); getErrorString();\n', "validar")
    r = ler_csv(os.path.join(PASTA, "MF_ST_res.csv"))
    tn = r["time"]
    degraus = np.array([-1.0, 1.0, 1.0, -1.0])  # h1 ... h4, os do MF_ST
    malhas = ("C1", "C2", "C3", "C4")

    i0 = np.searchsorted(tn, T_DEGRAU) - 1  # última amostra antes do degrau
    reg = {nm: r[f"P.{nm}"][i0] for nm in NIVEIS}
    print(f"\n7. Validação no MF_ST não linear (tanques partindo vazios, degraus de SP em t = {T_DEGRAU:.0f} s)")
    print("   regime antes do degrau:  " + "  ".join(f"{nm} = {v:.2f}" for nm, v in reg.items()))
    print("   aberturas:               " + "  ".join(f"{v} = {r[f'P.x{k + 1}'][i0]:.3f}" for k, v in enumerate(VALVULAS)))
    dif = max(abs(reg[nm] - h) for nm, h in zip(NIVEIS, h_op))
    print(f"   diferença para o ponto da linearização: {dif:.3f} cm")

    # grade uniforme depois do degrau (a saída do OpenModelica tem pontos extras nos eventos)
    tl = np.linspace(0, T_VALIDACAO, int(T_VALIDACAO) + 1)
    nl = {nm: np.interp(T_DEGRAU + tl, tn[i0 + 1:], r[nm][i0 + 1:]) for nm in
          [f"P.{h}" for h in NIVEIS] + [f"P.x{k}" for k in range(1, 5)] + [f"{c}.MV" for c in malhas]}
    _, yl, _ = signal.lsim(signal.StateSpace(Af, Bf, Cf, Df), np.tile(degraus, (len(tl), 1)), tl)
    fim = {nm: nl[f"P.{nm}"][-1] for nm in NIVEIS}
    print("   regime depois do degrau: " + "  ".join(f"{nm} = {v:.2f}" for nm, v in fim.items()))
    print("   aberturas:               " + "  ".join(f"{v} = {nl[f'P.x{k + 1}'][-1]:.3f}" for k, v in enumerate(VALVULAS)))
    baixos = [nm for nm in NIVEIS if min(reg[nm], fim[nm]) < 5]
    print("   todos os tanques >= 5 cm em regime, antes e depois" if not baixos else f"   ABAIXO DE 5 cm: {', '.join(baixos)}")
    xs = np.array([nl[f"P.x{k}"] for k in range(1, 5)])
    print(f"   aberturas durante o transitório: {xs.min():.3f} a {xs.max():.3f}")
    for k, nm in enumerate(NIVEIS[:4]):
        d = nl[f"P.{nm}"] - reg[nm]
        fora = np.nonzero(np.abs(d - degraus[k]) > 0.05)[0]
        pico = max(np.max(d * np.sign(degraus[k])) - 1, 0) * 100
        print(f"   {nm}: acomoda (±5 %) em {tl[fora[-1]] if len(fora) else 0:.0f} s, sobressinal {pico:.0f} %, "
              f"não linear x linear: {np.max(np.abs(d - yl[:, k])):.3f} cm")

    if sem_figura:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("   (sem matplotlib: figura não gerada)")
        return
    fig, ax = plt.subplots(2, 4, figsize=(15, 6), sharex=True)
    for k, nm in enumerate(NIVEIS):
        a = ax.flat[k]
        a.plot(tl, nl[f"P.{nm}"], label="não linear (MF_ST)")
        a.plot(tl, reg[nm] + yl[:, k], "--", label="linear")
        if min(nl[f"P.{nm}"]) < 6:
            a.axhline(5, color="0.6", lw=0.8, ls=":")
        a.set_title(nm)
        a.set_ylabel("cm")
    for a, grupo in zip(ax.flat[6:], (malhas[:2], malhas[2:])):
        for c in grupo:
            l = malhas.index(c)
            a.plot(tl, nl[f"{c}.MV"], label=f"{c}.MV")
            a.plot(tl, nl[f"{c}.MV"][0] + yl[:, 6 + l], "--", color=a.lines[-1].get_color())
        a.set_title("MV (abertura)")
        a.legend(fontsize=8)
    for a in ax[1]:
        a.set_xlabel(f"t - {T_DEGRAU:.0f} s")
    ax.flat[0].legend(fontsize=8)
    fig.suptitle("Degraus de SP (h1 -1, h2 +1, h3 +1, h4 -1 cm) com a sintonia SIMC: não linear (cheio) e linear (tracejado)")
    fig.tight_layout()
    saida = os.path.join(AQUI, "validacao.png")
    fig.savefig(saida, dpi=110)
    print(f"   figura: {os.path.relpath(saida, os.path.dirname(AQUI))}")

if __name__ == "__main__":
    main()
