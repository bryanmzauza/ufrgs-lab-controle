"""Servidor local do simulador de seis tanques.

O modelo é o próprio .mo: o OpenModelica compila o MF_ST uma vez e o executável
gerado roda em trechos curtos, um atrás do outro, sem tempo final. Cada trecho
começa do estado em que o anterior parou (-iif/-iit, a partir do resultado .mat
do trecho anterior) e usa os parâmetros do momento (-overrideFile). Assim a
simulação continua indefinidamente e responde aos ajustes da página. Nenhuma
equação do modelo é reescrita aqui nem na página.

Uso:
    python servidor.py                  # http://localhost:8000
    python servidor.py --porta 8080

Requer o OpenModelica (omc). O servidor procura, nesta ordem: a variável de
ambiente OMC, $OPENMODELICAHOME/bin, o PATH e C:/Program Files/OpenModelica*.
"""

import argparse
import glob
import hashlib
import json
import math
import os
import re
import shutil
import socket
import struct
import subprocess
import sys
import tempfile
import threading
import time
import xml.etree.ElementTree as ET
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

AQUI = os.path.dirname(os.path.abspath(__file__))
ARQUIVO_MO = os.path.normpath(os.path.join(AQUI, "..", "modelos", "modelo_6tanques_v11.mo"))
MODELO = "MF_ST"
PAGINA = os.path.join(AQUI, "simulador_seis_tanques.html")
PASTA_BUILD = os.path.join(tempfile.gettempdir(), "seis_tanques_om")
MAX_DURACAO = 3600      # s simulados por trecho
MAX_PONTOS = 20000      # pontos por trecho
TRECHOS_GUARDADOS = 40  # arquivos .mat mantidos para continuar a simulação

WINDOWS = sys.platform == "win32"


def achar_omc():
    candidatos = [os.environ.get("OMC")]
    home = os.environ.get("OPENMODELICAHOME")
    if home:
        candidatos.append(os.path.join(home, "bin", "omc.exe" if WINDOWS else "omc"))
    candidatos.append(shutil.which("omc"))
    if WINDOWS:
        candidatos += sorted(glob.glob(r"C:\Program Files\OpenModelica*\bin\omc.exe"), reverse=True)
    for c in candidatos:
        if c and os.path.isfile(c):
            return c
    return None


def ler_mat(caminho):
    """Lê o resultado .mat (MATLAB v4) do OpenModelica: {nome: [valores ao longo do tempo]}.

    Só as variáveis que variam no tempo (data_2); os parâmetros (data_1) ficam de fora.
    """
    with open(caminho, "rb") as f:
        buf = f.read()
    mats, pos = {}, 0
    tipos = {0: ("d", 8), 1: ("f", 4), 2: ("i", 4), 3: ("h", 2), 4: ("H", 2), 5: ("B", 1)}
    while pos < len(buf):
        tipo, linhas, colunas, _imag, tam_nome = struct.unpack_from("<5i", buf, pos)
        pos += 20
        nome = buf[pos:pos + tam_nome - 1].decode("latin-1")
        pos += tam_nome
        cod, tam = tipos[(tipo // 10) % 10]
        n = linhas * colunas
        mats[nome] = (linhas, colunas, struct.unpack_from(f"<{n}{cod}", buf, pos))
        pos += n * tam
    trans = bytes(mats["Aclass"][2][3::4]).decode("latin-1").strip().startswith("binTrans")
    if not trans:
        raise ErroModelo("Formato .mat inesperado (esperado binTrans).")
    nl, nv, chars = mats["name"]
    nomes = [bytes(chars[j * nl:(j + 1) * nl]).decode("latin-1").rstrip("\x00 ") for j in range(nv)]
    info = mats["dataInfo"][2]
    n2, nt, d2 = mats["data_2"]
    saida = {}
    for j, nome in enumerate(nomes):
        qual, idx = info[4 * j], info[4 * j + 1]
        if qual not in (0, 2) or nome.startswith("$"):
            continue
        linha, sinal = abs(idx) - 1, (-1.0 if idx < 0 else 1.0)
        saida[nome] = [sinal * d2[c * n2 + linha] for c in range(nt)]
    return saida


class ErroModelo(Exception):
    """Falha ao compilar ou simular; a mensagem vai para a página."""


class TrechoPerdido(Exception):
    """O trecho de onde a página quer continuar não existe mais (modelo recompilado ou servidor reiniciado)."""


class Modelo:
    def __init__(self, omc):
        self.omc = omc
        self.lock = threading.Lock()
        self.hash = None
        self.pasta = None
        self.info = None
        self.trechos = {}   # id -> arquivo .mat, na ordem do uso mais recente
        self.proximo = 1

    # ---------- compilação ----------
    def garantir(self):
        """Compila o .mo se ele mudou desde a última compilação."""
        with open(ARQUIVO_MO, "rb") as f:
            fonte = f.read()
        h = hashlib.sha256(fonte).hexdigest()[:16]
        if h == self.hash:
            return
        pasta = os.path.join(PASTA_BUILD, h)
        if not os.path.isfile(os.path.join(pasta, MODELO + "_init.xml")) or not self._executavel(pasta):
            self._compilar(pasta)
        self.info = self._ler_parametros(pasta, fonte.decode("utf-8", errors="replace"), h)
        self.hash, self.pasta = h, pasta
        # os trechos do modelo anterior não servem para continuar o novo
        for arq in self.trechos.values():
            if os.path.exists(arq):
                os.remove(arq)
        self.trechos = {}

    def _executavel(self, pasta):
        # No Windows o omc gera o .bat (que põe as DLLs do OpenModelica no PATH) antes de compilar o C,
        # então o .bat sozinho não prova que a compilação deu certo: o .exe precisa existir.
        exe = os.path.join(pasta, MODELO + (".exe" if WINDOWS else ""))
        if not os.path.isfile(exe):
            return None
        bat = os.path.join(pasta, MODELO + ".bat")
        return bat if WINDOWS and os.path.isfile(bat) else exe

    def _compilar(self, pasta):
        shutil.rmtree(pasta, ignore_errors=True)
        os.makedirs(pasta)
        caminho_mo = ARQUIVO_MO.replace("\\", "/")
        script = (
            f'loadFile("{caminho_mo}"); getErrorString();\n'
            f'buildModel({MODELO}, outputFormat="mat"); getErrorString();\n'
        )
        with open(os.path.join(pasta, "build.mos"), "w", encoding="utf-8") as f:
            f.write(script)
        t0 = time.time()
        r = subprocess.run([self.omc, "build.mos"], cwd=pasta, capture_output=True, text=True, timeout=600)
        saida = (r.stdout + r.stderr).strip()
        if not self._executavel(pasta):
            shutil.rmtree(pasta, ignore_errors=True)
            raise ErroModelo("O OpenModelica não conseguiu compilar o modelo.\n\n" + saida)
        with open(os.path.join(pasta, "build.log"), "w", encoding="utf-8") as f:
            f.write(f"{time.time() - t0:.1f} s\n{saida}\n")

    def _ler_parametros(self, pasta, fonte, h):
        linhas = fonte.splitlines()
        raiz = ET.parse(os.path.join(pasta, MODELO + "_init.xml")).getroot()
        params = []
        for sv in raiz.iter("ScalarVariable"):
            real = sv.find("Real")
            if sv.get("causality") != "parameter" or sv.get("isValueChangeable") != "true":
                continue
            if real is None or real.get("start") is None:
                continue
            nome = sv.get("name")
            params.append({
                "nome": nome,
                "valor": float(real.get("start")),
                "linha": int(sv.get("startLine") or 0),
                "comentario": self._comentario(linhas, nome.split(".")[-1], sv),
            })
        versao = subprocess.run([self.omc, "--version"], capture_output=True, text=True).stdout.strip()
        return {
            "modelo": MODELO,
            "arquivo": os.path.relpath(ARQUIVO_MO, os.path.dirname(AQUI)).replace("\\", "/"),
            "hash": h,
            "openmodelica": versao,
            "parametros": params,
        }

    @staticmethod
    def _comentario(linhas, curto, sv):
        """Comentário // da linha em que o parâmetro é declarado no .mo."""
        a, b = int(sv.get("startLine") or 0), int(sv.get("endLine") or 0)
        padrao = re.compile(r"(?<![\w.])" + re.escape(curto) + r"\s*(=|\(|,|;)")
        for n in range(a, b + 1):
            if 0 < n <= len(linhas) and padrao.search(linhas[n - 1].split("//")[0]):
                partes = linhas[n - 1].split("//", 1)
                return partes[1].strip() if len(partes) > 1 else ""
        return ""

    # ---------- simulação em trechos ----------
    def trecho(self, alterados, de, duracao, passo):
        """Simula [t0, t0 + duracao]. Sem "de", parte de t = 0 com a inicialização do .mo;
        com "de" = {"trecho": id, "t": t0}, continua do estado daquele trecho no instante t0."""
        with self.lock:
            hash_antes = self.hash
            self.garantir()
            if de and self.hash != hash_antes and hash_antes is not None:
                raise TrechoPerdido("O .mo mudou e foi recompilado.")
            validos = {p["nome"] for p in self.info["parametros"]}
            for nome, v in alterados.items():
                if nome not in validos:
                    raise ErroModelo(f"Parâmetro desconhecido no {MODELO}: {nome}")
                if not isinstance(v, (int, float)) or not math.isfinite(v):
                    raise ErroModelo(f"Valor inválido para {nome}: {v!r}")
            if not (0 < duracao <= MAX_DURACAO and passo > 0):
                raise ErroModelo("Duração do trecho ou passo de saída inválidos.")
            passo = max(passo, duracao / MAX_PONTOS)
            t0, extra = 0.0, []
            if de:
                arq = self.trechos.pop(int(de["trecho"]), None)
                if not arq or not os.path.exists(arq):
                    raise TrechoPerdido("O trecho anterior não está mais no servidor.")
                self.trechos[int(de["trecho"])] = arq  # volta ao fim da fila: usado por último, descartado por último
                t0 = float(de["t"])
                # Sem -iim=none de propósito: a inicialização precisa rodar para recalcular o que o
                # OpenModelica deriva só de parâmetros (ex.: P.rot_bomba_inicial = rot_bomba_inicial_bias).
                # Pulando-a, esses valores vinham do trecho anterior e o parâmetro novo não tinha efeito.
                extra = [f"-iif={os.path.basename(arq)}", f"-iit={t0!r}"]
            pasta = self.pasta
            with open(os.path.join(pasta, "override.txt"), "w", encoding="utf-8") as f:
                for nome, v in alterados.items():
                    f.write(f"{nome}={float(v)!r}\n")
            # startTime, stopTime e stepSize não aceitam -override: vão numa cópia do XML de inicialização
            xml = open(os.path.join(pasta, MODELO + "_init.xml"), encoding="utf-8").read()
            for campo, valor in (("startTime", t0), ("stopTime", t0 + duracao), ("stepSize", passo), ("outputFormat", "mat")):
                xml = re.sub(rf'({campo}\s*=\s*")[^"]*"', lambda m: f'{m.group(1)}{valor}"', xml, count=1)
            with open(os.path.join(pasta, "run_init.xml"), "w", encoding="utf-8") as f:
                f.write(xml)
            ident = self.proximo
            self.proximo += 1
            saida = f"trecho_{ident}.mat"
            args = [self._executavel(pasta), "-f=run_init.xml", f"-r={saida}"] + extra
            if alterados:  # o runtime recusa um overrideFile vazio
                args.append("-overrideFile=override.txt")
            t_ini = time.time()
            r = subprocess.run(args, cwd=pasta, capture_output=True, text=True, timeout=300)
            dt = time.time() - t_ini
            log = (r.stdout + r.stderr).strip()
            caminho = os.path.join(pasta, saida)
            if r.returncode != 0 or not os.path.exists(caminho):
                raise ErroModelo("A simulação falhou.\n\n" + log)
            dados = ler_mat(caminho)
            self.trechos[ident] = caminho
            for velho in list(self.trechos)[:-TRECHOS_GUARDADOS]:  # ordem de uso: descarta os menos recentes
                if os.path.exists(self.trechos[velho]):
                    os.remove(self.trechos[velho])
                del self.trechos[velho]
            avisos = [l.split("|", 2)[-1].strip() for l in log.splitlines() if "| warning |" in l or "| error" in l]
            t = dados.pop("time")
            return {
                "hash": self.hash,
                "trecho": ident,
                "segundos": round(dt, 3),
                "avisos": avisos,
                "t": [round(x, 6) for x in t],
                "vars": {n: [float(f"{x:.6g}") for x in c] for n, c in dados.items()},
            }

    def descrever(self):
        with self.lock:
            self.garantir()
            return self.info


class Handler(BaseHTTPRequestHandler):
    modelo = None

    def log_message(self, fmt, *args):
        pass  # um pedido por trecho: o log do terminal ficaria ilegível

    def _json(self, codigo, obj):
        corpo = json.dumps(obj, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        self.send_response(codigo)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(corpo)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(corpo)

    def do_GET(self):
        caminho = self.path.split("?")[0]
        if caminho in ("/", "/simulador_seis_tanques.html"):
            with open(PAGINA, "rb") as f:
                corpo = f.read()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(corpo)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(corpo)
        elif caminho == "/api/modelo":
            self._executar(self.modelo.descrever)
        else:
            self.send_error(404)

    def do_POST(self):
        if self.path.split("?")[0] != "/api/trecho":
            self.send_error(404)
            return
        try:
            n = int(self.headers.get("Content-Length") or 0)
            pedido = json.loads(self.rfile.read(n) or b"{}")
            alterados = pedido.get("parametros") or {}
            de = pedido.get("de")
            duracao = float(pedido["duracao"])
            passo = float(pedido["passo"])
        except (KeyError, ValueError, TypeError, json.JSONDecodeError) as e:
            self._json(400, {"erro": f"Pedido inválido: {e}"})
            return
        self._executar(lambda: self.modelo.trecho(alterados, de, duracao, passo))

    def _executar(self, fn):
        try:
            self._json(200, fn())
        except TrechoPerdido as e:
            self._json(409, {"erro": str(e), "reiniciar": True})
        except ErroModelo as e:
            self._json(422, {"erro": str(e)})
        except subprocess.TimeoutExpired:
            self._json(504, {"erro": "O OpenModelica excedeu o tempo limite."})
        except Exception as e:  # noqa: BLE001 - qualquer outra falha vira mensagem na página
            self._json(500, {"erro": f"{type(e).__name__}: {e}"})


def main():
    ap = argparse.ArgumentParser(description="Servidor do simulador de seis tanques (OpenModelica).")
    ap.add_argument("--host", default=None, help="padrão: só a máquina local (127.0.0.1 e ::1)")
    ap.add_argument("--porta", type=int, default=8000)
    a = ap.parse_args()

    omc = achar_omc()
    if not omc:
        sys.exit("OpenModelica (omc) não encontrado. Instale o OpenModelica ou defina a variável OMC com o caminho do omc.")
    Handler.modelo = Modelo(omc)
    print(f"omc: {omc}")
    print(f"modelo: {ARQUIVO_MO} ({MODELO})")
    print("compilando...", flush=True)
    try:
        info = Handler.modelo.descrever()
        print(f"pronto: {len(info['parametros'])} parâmetros, {info['openmodelica']}")
    except ErroModelo as e:
        print(e, file=sys.stderr)
        print("O servidor sobe mesmo assim; corrija o .mo e recarregue a página.", file=sys.stderr)
    srv = ThreadingHTTPServer((a.host or "127.0.0.1", a.porta), Handler)
    if a.host is None and socket.has_ipv6:
        # No Windows, "localhost" tenta ::1 antes de 127.0.0.1; sem ninguém escutando em ::1,
        # cada pedido esperaria ~2 s pelo fallback. Escutar nos dois evita essa espera.
        class ServidorV6(ThreadingHTTPServer):
            address_family = socket.AF_INET6
        try:
            v6 = ServidorV6(("::1", a.porta), Handler)
            threading.Thread(target=v6.serve_forever, daemon=True).start()
        except OSError:
            pass
    print(f"abra http://localhost:{a.porta}", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
