# Bancada de seis tanques: Laboratório de Controle

Simulador no navegador do sistema de seis tanques com quatro malhas de nível, desenvolvido para a disciplina de Laboratório de Controle (mestrado em Engenharia Química). A página roda o modelo em Modelica no OpenModelica.

## Estrutura

```
modelos/
  modelo_6tanques_v11.mo        modelo (SixTanks, PID_ISA, MF_ST)
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
- O servidor compila o `MF_ST` de `modelos/modelo_6tanques_v11.mo` com o `omc`, uma vez só.
- O executável roda em trechos curtos, sem tempo final. Cada trecho continua do estado em que o anterior parou (`-iif` do OpenModelica).
- Os sliders mudam parâmetros do `.mo` com a simulação rodando: o próximo trecho já usa o valor novo.
- Quando o `.mo` muda, o servidor recompila sozinho e a simulação reinicia.
- A página não tem nenhuma equação do modelo: só envia os parâmetros e desenha o resultado.

A página segue o padrão de tela de SDCD (ISA-101): fundo cinza, cor só para alarme e barra de alarmes sempre visível no rodapé. Ela tem cinco telas, pensadas para caber numa tela Full HD sem rolagem:
- **Visão geral:** sinótico, tendências e um faceplate por malha (LIC-01 a LIC-04) com barras de PV/SP/MV. O SP é digitado (Enter aplica) ou ajustado em ▲/▼. Clicar num tanque chama o controlador dele;
- **Parâmetros:** sliders para os parâmetros do `MF_ST` (controladores, processo, bombas) e campos para os demais;
- **Análise:** vazões, pressões, estado atual e todas as variáveis do resultado;
- **Alarmes:** alarmes ativos com reconhecimento e diário de eventos (alarmes, reconhecimentos e ações do operador);
- **Cenários:** 5 conjuntos de parâmetros prontos para apresentar.

Os alarmes são da página, não do `.mo`: comparam o resultado do OpenModelica com limites de operação (nível alto, nível abaixo de 5 cm em regime, desvio PV − SP, válvula saturada e falha de comunicação com o OpenModelica).

O mesmo `.mo` também abre no OMEdit ou no Dymola: simule `MF_ST`.
