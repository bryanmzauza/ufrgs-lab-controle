# Bancada de seis tanques: Laboratório de Controle

Simulador no navegador do sistema de seis tanques com quatro malhas de nível, desenvolvido para a disciplina de Laboratório de Controle (mestrado em Engenharia Química). A página roda o modelo em Modelica no OpenModelica.

## Estrutura

```
modelos/
  modelo_6tanques_v12.mo        modelo (SixTanks, PID_ISA, MF_ST e SixTanks_MA)
  historico/                    versões anteriores do modelo
linearizacao/
  analise.py                    linearização no ponto de operação, RGA e sintonia SIMC
simulador/
  servidor.py                   compila o .mo com o OpenModelica e simula a pedido da página
  simulador_seis_tanques.html   página servida pelo servidor.py (gráficos, diagrama, parâmetros)
```

## Como usar

Requer o [OpenModelica](https://openmodelica.org) e o Python 3.

```
cd simulador
python servidor.py
```

Depois, abra `http://localhost:8000`.

Como funciona:
- O servidor compila o `MF_ST` de `modelos/modelo_6tanques_v12.mo` com o `omc`, uma vez só.
- O executável roda em trechos curtos, sem tempo final. Cada trecho continua do estado em que o anterior parou (`-iif` do OpenModelica).
- Os sliders mudam parâmetros do `.mo` com a simulação rodando: o próximo trecho já usa o valor novo.
- Quando o `.mo` muda, o servidor recompila sozinho e a simulação reinicia.
- A página não tem nenhuma equação do modelo: só envia os parâmetros e desenha o resultado.

A página segue o padrão de tela de SDCD (ISA-101): fundo cinza, cor só para alarme e barra de alarmes sempre visível no rodapé. Ela tem cinco telas, pensadas para caber numa tela Full HD sem rolagem:
- **Visão geral:** sinótico, tendências e um faceplate por malha (LIC-01 a LIC-04) com barras de PV/SP/MV. O SP é digitado (Enter aplica) ou ajustado em ▲/▼. Clicar num tanque chama o controlador dele;
- **Parâmetros:** sliders para os parâmetros do `MF_ST` (controladores, processo, bombas) e campos para os demais;
- **Análise:** vazões, pressões, estado atual e todas as variáveis do resultado;
- **Alarmes:** alarmes ativos com reconhecimento e diário de eventos (alarmes, reconhecimentos e ações do operador);
- **Cenários:** 6 conjuntos de parâmetros prontos para apresentar: o `.mo` como está e cinco variações da v11.

Os alarmes são da página, não do `.mo`: comparam o resultado do OpenModelica com limites de operação (nível alto, nível abaixo de 5 cm em regime, desvio PV − SP, válvula saturada e falha de comunicação com o OpenModelica).

O mesmo `.mo` também abre no OMEdit ou no Dymola: simule `MF_ST`.

## Linearização e sintonia

```
python linearizacao/analise.py
```

Requer numpy e scipy (e matplotlib para a figura). O script:
1. lineariza a `SixTanks_MA` com o `linearize()` do OpenModelica. A `SixTanks_MA` é a `SixTanks` em malha aberta, que se inicializa no equilíbrio do ponto de operação do `MF_ST`;
2. mostra A, B, as funções de transferência, o RGA e os zeros de transmissão;
3. calcula a sintonia SIMC (PI): primeiro C3/C4 (T3 e T4 são integradores), depois C1/C2 com C3/C4 fechados;
4. valida no `MF_ST` não linear: aplica os degraus de SP a partir do regime, compara com o modelo linear e salva `linearizacao/validacao.png`.

Resultados no ponto de operação da v12 (B-01 a 100 %, r1 = r2 = 0,15, SP de h1/h2 em 9 cm e de h3/h4 em 12 cm):
- os seis tanques ficam acima de 5 cm (h5 = h6 = 6,9 cm) e as válvulas ficam em 57 % (FV-01/02) e 24 % (FV-03/04);
- RGA com C3/C4 fechados: λ11 = −0,03 e λ12 = 1,03, então o pareamento é cruzado (C1 → FV-02, C2 → FV-01). Com r1 + r2 < 1 há um zero de transmissão no semiplano direito, em +0,18 1/s;
- sintonia: C1/C2 com Kp = 0,062 e Ti = 33 s; C3/C4 com Kp = −0,051 e Ti = 40 s (Td = 0);
- com degraus de 1 cm, as malhas acomodam em 2 a 3 min, com sobressinal de até 8 %. O modelo não linear e o linear diferem em até 0,12 cm.

As equações da planta (`SixTanks`) são as mesmas da v11. O cabeçalho do `.mo` explica cada mudança.
