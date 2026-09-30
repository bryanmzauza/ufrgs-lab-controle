// =====================================================================
// Modelo seis tanques -- versão 12
// ---------------------------------------------------------------------
// Mudanças em relação à v11.
// As linhas alteradas estão marcadas com "v12:". As marcas "v11:" são da versão anterior.
//
//  A física não mudou: SixTanks, funções auxiliares e PID_ISA são idênticos
//  aos da v11. Mudaram o ponto de operação, o pareamento e a sintonia do MF_ST,
//  e entrou a SixTanks_MA (planta em malha aberta para linearização).
//
//  1. Ponto de operação viável (MF_ST). Na v11, com os parâmetros padrão,
//     FV-01 e FV-02 ficavam saturadas em 100% e o regime era h1 = 3,6 cm e
//     h2 = 1,8 cm, com SP de 11 e 13 cm. Dois motivos:
//     - B-01 a 40% entrega 15,9 cm³/s por ramo, e T1/T2 a 12 cm pedem
//       CD*12^0,6 = 39 cm³/s cada. Mesmo a 100% a B-01 dá no máximo 39,8 cm³/s
//       por ramo (divisão de 50%);
//     - com r1 + r2 = 1, T3 recebe sempre r1*(F1 + F2) e T4 (1 - r1)*(F1 + F2):
//       o ganho estático de (h1, h2) em relação a (x1, x2) é singular e
//       h1/h2 fica fixo (~1,97 com r1 = 0,6), qualquer que seja a sintonia.
//     O requisito de T5/T6 >= 5 cm em regime pede (1 - r)*F >= CD*5^0,6 = 23,3 cm³/s,
//     então r1, r2 <= ~0,3. Ficou: rot_bomba_inicial_bias = 1, r1 = r2 = 0,15,
//     SP de h1/h2 = 9 cm (T3/T4 continuam em 12). No regime, h5 = h6 = 6,9 cm e as
//     válvulas ficam em 57% (FV-01/02) e 24% (FV-03/04). Depois dos degraus de SP,
//     h5 = 5,8 cm e h6 = 8,0 cm, e as válvulas ficam entre 19% e 70% durante o transitório.
//
//  2. Pareamento cruzado: C1 (h1) -> FV-02, C2 (h2) -> FV-01. Com r < 0,5 a maior
//     parte de F2 chega a T1 por T5 -> T3 -> T1. RGA com C3/C4 fechados:
//     λ11 = -0,03 (pareamento da v11) e λ12 = 1,03. O parâmetro pareamento_cruzado = 0
//     volta ao pareamento da v11. r1 + r2 < 1 também dá um zero de transmissão no
//     semiplano direito (+0,18 1/s), longe da banda das malhas de nível.
//
//  3. Sintonia SIMC (PI, Td = 0) a partir do modelo linearizado:
//     - C3/C4: T3/T4 são integradores (a vazão de B-03/B-04 não depende do nível):
//       G = -1,95/s, τc = 10 s -> Kp = -0,051, Ti = 40 s;
//     - C1/C2: com C3/C4 fechados, G ≈ 12,0 e^(-11,5s)/(22,3s + 1)² -> regra da
//       metade e τc = θ -> Kp = 0,062, Ti = 33 s.
//     Validação no não linear (degraus de 1 cm a partir do regime): acomoda em
//     2 a 3 min, sobressinal <= 8% e diferença não linear x linear <= 0,12 cm.
//     Na v11 as quatro malhas usavam Kp = ±0,05, Ti = 120 s, Td = 1 s.
//     TimeYset passou de 10 s para 300 s: com 10 s o degrau caía no meio do
//     enchimento; em 300 s os níveis já estão no SP e dá para ver a resposta servo.
//
//  4. SixTanks_MA (fim do arquivo): planta em malha aberta que se inicializa
//     no equilíbrio do ponto de operação, para linearize(). O cálculo está em
//     linearizacao/analise.py.
//
//  Observações que continuam valendo (não são correções):
//  - PID_ISA: Id0 = Yset_bias supõe o nível já no setpoint. Com Td = 0 a parcela
//    derivativa fica desligada, então o pico Ud(0) na partida da v11 não aparece mais.
//  - CD = 0.20*sqrt(2*981) vem de Torricelli (expoente 0,5), mas é usado com
//    expoente 0,6: na prática é um coeficiente empírico, a ajustar com a bancada.
//    Se mudar, o ponto de operação e a sintonia precisam ser refeitos (analise.py).
//  - ΔP não tem unidade declarada (K_bomba = 100, CV = 15).
// =====================================================================


// ---------- funções auxiliares (v11) ----------

function vazao_descarga "Descarga por gravidade q = CD*h^p, que vai a zero com o tanque vazio"
  input Real h "nível (cm)";
  input Real CD "coeficiente de descarga";
  input Real p "expoente de descarga";
  input Real h_reg = 1e-3 "abaixo deste nível a curva vira uma reta até zero (h^p tem derivada infinita em h = 0)";
  output Real q "vazão (cm³/s)";
algorithm
  q := if noEvent(h >= h_reg) then CD*h^p
       elseif noEvent(h > 0) then CD*h_reg^(p - 1)*h
       else 0;
end vazao_descarga;

function raiz_reg "Raiz quadrada que vai a zero: sqrt(x) para x >= eps, reta até zero abaixo disso"
  input Real x;
  input Real eps = 1e-6;
  output Real y;
algorithm
  y := if noEvent(x >= eps) then sqrt(x)
       elseif noEvent(x > 0) then x/sqrt(eps)
       else 0;
end raiz_reg;

function fator_succao "Fração da vazão que a bomba consegue puxar, pelo nível do tanque de sucção"
  input Real h "nível do tanque de sucção (cm)";
  input Real h_min "a partir deste nível a bomba tem sucção plena (cm)";
  output Real f "0 com o tanque vazio, 1 com h >= h_min";
algorithm
  f := if noEvent(h >= h_min) then 1
       elseif noEvent(h > 0) then h/h_min
       else 0;
end fator_succao;


model PID_ISA


  // Parametros do Controlador
  parameter Real Kp=1, Ti=15, Td=1, b=1, c=0, N=10, FTt = 2;
  parameter Real U_bias=0.5, Umax=1, Umin=0, Yset_bias=0; // edit: Umax e Umin e Ubias = 1, 0, 0.5 para refletir uma abertura de válvula de controle
  
  // Condicoes Iniciais
  parameter Real Ui0 = Kp * Yset_bias * (1-b), 
                 Id0 = (Yset_bias - c*Yset_bias);
  
  // Variaveis do Controlador
  Real Ui (start=Ui0), Id (start= Id0);
  Real Uop (start=U_bias),  Ud(start=0), Up (start= 0);
 
  // Inputs e Outputs
  input  Real SP (start = Yset_bias), PV (start= Yset_bias);
  output Real MV (start = U_bias);

equation
     // Parcela Integral -- Ui
     Ti*der(Ui) = Kp*(SP-PV)+FTt*(MV - Uop);
     // Parcela Derivativa  -- Ud
     der(Id) = (if Td > 0.0  then -N*(Id+c*SP-PV)/Td else 0);
     Ud = if Td > 0.0  then Kp*N*(Id+c*SP-PV) else 0;
     // Parcela Proporcional -- Up
     Up = Kp*(b*SP-PV);      
     // Acao total de controle calculada pelo controlador -- Uop -- OP (OutPut)
     Uop = U_bias + Up + Ui + Ud;
     // Acao de controle efetivamente aplicada na planta -- MV -- manipulada
     MV = smooth(0,if Uop > Umax then Umax else if Uop < Umin then Umin else Uop);
end PID_ISA;

// ==========================================
// ----- Modelo seis tanques          ------
// ==========================================
model SixTanks

    // parâmetros das descargas dos tanques
    // v11: CD3 e CD4 removidos (T3 e T4 não descarregam por gravidade; quem os esvazia são B-03 e B-04)
    parameter Real CD1=0.20*sqrt(2*981), CD2=0.20*sqrt(2*981),
                   CD5=0.20*sqrt(2*981), CD6=0.20*sqrt(2*981);

    parameter Real A1 = 63, A2 = 63, A3 = 63, A4 =63, A5=63, A6=63;  // cm²

    parameter Real F_max_ramo = 83.0; // vazao máxima que a bomba consegue fornecer para cada um dos lados (entrada de água no sistema) (5L/min -> cm³/s)
    parameter Real F_max_bomba_saida3 = 83.0; // vazao máxima que a bomba consegue puxar (saída de água tanque 3)   (5L/min -> cm³/s)  -- v11: era Real, virou parameter
    parameter Real F_max_bomba_saida4 = 83.0; // vazao máxima que a bomba consegue puxar (saída de água tanque 4)   (5L/min -> cm³/s)  -- v11: era Real, virou parameter


    // potencia de descarga dos tanques. não será 0.5 pq não e totalmente turbulento
    // tb não será 1 pq não é totalmente laminar. 0.5 <= potencia <= 1.0
    // v11: pot_descarga3 e pot_descarga4 removidos (sem uso, como CD3 e CD4)
    parameter Real pot_descarga1 = 0.6, pot_descarga2 = 0.6,
                   pot_descarga5 = 0.6, pot_descarga6 = 0.6;

    parameter Real h_succao_min = 0.5; // v11: nível (cm) abaixo do qual B-03/B-04 começam a perder sucção; com o tanque vazio a vazão é zero


    // Variáveis do sistema
    Real F_bomba_inicial(start=10.0), F1(start=5.0), F2(start=5.0);
    Real F1_completo(start=10.0), F2_completo(start=10.0);
    Real F3(start=5.0), F4(start=5.0);
    Real F1_direto, F1_indireto, F2_direto, F2_indireto;
    Real F3_valvula_controle(start=4.0), F3_refluxo(start=1.0);
    Real F4_valvula_controle(start=4.0), F4_refluxo(start=1.0);
    Real q1, q2, q5, q6;                   // v11: descargas por gravidade (cm³/s), antes só dentro das equações
    Real fator_succao_3, fator_succao_4;   // v11: 0 a 1, sucção disponível em B-03/B-04

    Real funcao_obturador_valvula_1, funcao_obturador_valvula_2;
    Real funcao_obturador_valvula_3, funcao_obturador_valvula_4;

    Real delta_P_max_inicial(start=50.0);  // v11: substitui delta_P_max_esq e delta_P_max_dir, que tinham a mesma equação
    Real delta_P_bomba_lado_esq(start=40.0), delta_P_bomba_lado_dir(start=40.0);
    Real delta_P_max_3(start=50.0);
    Real delta_P_max_4(start=50.0);
    Real delta_P_bomba_3(start=40.0), delta_P_bomba_4(start=40.0);
    Real F_max_atual, F_max_atual_bomba3, F_max_atual_bomba4;


    Real h1(start=0, fixed=true), h2(start=0, fixed=true),     // nível inicial em cm (0 = tanque vazio)
         h3(start=0, fixed=true), h4(start=0, fixed=true),     // nível inicial em cm (0 = tanque vazio)
         h5(start=0, fixed=true), h6(start=0, fixed=true);     // nível inicial em cm (0 = tanque vazio)


    // criando variáveis para atuar -- abertura das válvulas
    input Real x1; // Abertura da FV-01, após a bomba B-01 (Perna Esquerda, lado ímpar)
    input Real x2; // Abertura da FV-02, após a bomba B-01 (Perna Direita, lado par)
    input Real x3; // Abertura da FV-03, saída do T3
    input Real x4; // Abertura da FV-04, saída do T4

    input   Real r1, r2;    // divisão das correntes de entrada (valvula de 3 vias)

    input   Real rot_bomba_inicial, rot_bomba3, rot_bomba4; // Rotações (0 a 1)

    input   Real divisao_bomba_inicial; // Multiplexação temporal (0 a 1) -- % do tempo que vai para cada lado

    output  Real h1out, h2out, h3out, h4out, h5out, h6out;  // níveis do tanque. variáveis diferentes para caso queira inserir um atraso de medida


    parameter Real K_bomba_inicial=100; // parametros da equacao do delta P maximo da bomba
    parameter Real K_bomba_3=100; // parametros da equacao do delta P maximo da bomba
    parameter Real K_bomba_4=100; // parametros da equacao do delta P maximo da bomba


    parameter Real CV_valvula1=15, CV_valvula2=15, CV_valvula3=15, CV_valvula4=15;
    parameter Real abertura_valv_refluxo_3=0.2, abertura_valv_refluxo_4=0.2;
    parameter Real CV_valvula3_refluxo=5, CV_valvula4_refluxo=5;


equation

    F_max_atual = rot_bomba_inicial * F_max_ramo; // considerando, inicialmente, que a vazão máxima varia linearmente com a rotação da bomba

    funcao_obturador_valvula_1 = max(0.0, x1); // obturador linear por enquanto
    funcao_obturador_valvula_2 = max(0.0, x2); // obturador linear por enquanto

    // BOMBA B-01 - PUXA ÁGUA DO RESERVATÓRIO

    // ideia será calcular qual a vazão em cada ramo se fosse só por aquele ramo e ponderar por quanto tempo ele vai para cada ramo
    // variável "completo" é caso fosse tudo para aquele lado

    delta_P_max_inicial = K_bomba_inicial * rot_bomba_inicial^2;  // Delta P max é função da rotação da bomba (v11: um só para os dois lados)

    // lado esquerdo (FV-01)
    delta_P_bomba_lado_esq = delta_P_max_inicial * sqrt(max(1e-6, 1.0 - (F1_completo / max(1e-5, F_max_atual))^2));  // eq da bomba
    F1_completo = funcao_obturador_valvula_1 * CV_valvula1 * raiz_reg(delta_P_bomba_lado_esq);  // v11: raiz_reg no lugar de sqrt(max(1e-6, .)): zero com a bomba parada


    // lado direito (FV-02)
    delta_P_bomba_lado_dir = delta_P_max_inicial * sqrt(max(1e-6, 1.0 - (F2_completo / max(1e-5, F_max_atual))^2));  // eq da bomba
    F2_completo = funcao_obturador_valvula_2 * CV_valvula2 * raiz_reg(delta_P_bomba_lado_dir);  // v11: raiz_reg

    F1 = F1_completo * divisao_bomba_inicial;       // ponderação pelo tempo que joga água para a esq
    F2 = F2_completo * (1 - divisao_bomba_inicial); // ponderação pelo tempo que joga água para a dir

    F_bomba_inicial = F1 + F2; // conservação de massa


    // --- Divisão Válvulas de 3 Vias --- r é o quanto vai para o tanque direto
    F1_direto = r1 * F1;          // quanto vai de F1 para o tanque 3
    F1_indireto = (1 - r1) * F1;  // quanto vai de F1 para o tanque 6

    F2_direto = r2 * F2;          // quanto vai de F2 para o tanque 4
    F2_indireto = (1 - r2) * F2;  // quanto vai de F2 para o tanque 5


    // descargas por gravidade -- v11: vazao_descarga vai a zero com o tanque vazio (v10: CD*max(h, 1e-4)^p)
    q1 = vazao_descarga(h1, CD1, pot_descarga1);
    q2 = vazao_descarga(h2, CD2, pot_descarga2);
    q5 = vazao_descarga(h5, CD5, pot_descarga5);
    q6 = vazao_descarga(h6, CD6, pot_descarga6);

    A5 * der ( h5 ) = F2_indireto - q5;  // descarga por gravidade
    A6 * der ( h6 ) = F1_indireto - q6;  // descarga por gravidade



    // tanques intermediarios - tem a valvula intermediaria e bomba

    // BOMBAS B-03 E B-04 -- INTERMEDIARIAS

    F_max_atual_bomba3 = rot_bomba3 * F_max_bomba_saida3; // considerando, inicialmente, que a vazão máxima varia linearmente com a rotação da bomba
    F_max_atual_bomba4 = rot_bomba4 * F_max_bomba_saida4; // considerando, inicialmente, que a vazão máxima varia linearmente com a rotação da bomba

    funcao_obturador_valvula_3 = max(0.0, x3); // obturador linear por enquanto
    funcao_obturador_valvula_4 = max(0.0, x4); // obturador linear por enquanto


    delta_P_max_3 = K_bomba_3 * rot_bomba3^2;    // Delta P max é função da rotação da bomba
    delta_P_bomba_3 = delta_P_max_3 * sqrt(max(1e-6, 1.0 - (F3 / max(1e-5, F_max_atual_bomba3))^2));             // curva da bomba

    delta_P_max_4 = K_bomba_4 * rot_bomba4^2;    // Delta P max é função da rotação da bomba
    delta_P_bomba_4 = delta_P_max_4 * sqrt(max(1e-6, 1.0 - (F4 / max(1e-5, F_max_atual_bomba4))^2));             // curva da bomba

    // v11: sucção de B-03/B-04 pelo nível de T3/T4 (na v10 a bomba puxava água de tanque vazio)
    fator_succao_3 = fator_succao(h3, h_succao_min);
    fator_succao_4 = fator_succao(h4, h_succao_min);

    // Equacoes de saída

    F3 = F3_refluxo + F3_valvula_controle;
    F3_valvula_controle = fator_succao_3 * funcao_obturador_valvula_3 * CV_valvula3 * raiz_reg(delta_P_bomba_3);  // v11: fator_succao e raiz_reg
    F3_refluxo = fator_succao_3 * CV_valvula3_refluxo * abertura_valv_refluxo_3 * raiz_reg(delta_P_bomba_3);      // v11: fator_succao e raiz_reg


    F4 = F4_refluxo + F4_valvula_controle;
    F4_valvula_controle = fator_succao_4 * funcao_obturador_valvula_4 * CV_valvula4 * raiz_reg(delta_P_bomba_4);  // v11: fator_succao e raiz_reg
    F4_refluxo = fator_succao_4 * CV_valvula4_refluxo * abertura_valv_refluxo_4 * raiz_reg(delta_P_bomba_4);      // v11: fator_succao e raiz_reg


    // BALANÇO TANQUES INTERMEDIARIOS E SAÍDA

    A4 * der ( h4 ) = q6 + F2_direto - F4_valvula_controle; // entra o que sai do tanque 6 e o desvio direto de F2; sai o que a FV-04 deixa passar (o refluxo sai e volta ao próprio tanque)
    A3 * der ( h3 ) = q5 + F1_direto - F3_valvula_controle; // entra o que sai do tanque 5 e o desvio direto de F1; sai o que a FV-03 deixa passar (o refluxo sai e volta ao próprio tanque)

    A2 * der ( h2 ) = F4_valvula_controle - q2;
    A1 * der ( h1 ) = F3_valvula_controle - q1;

    h1out = h1; h2out = h2; h3out = h3;
    h4out = h4; h5out = h5; h6out = h6;

end SixTanks;
// v11: removido o comentário com o dicionário S1 (valores antigos que não correspondiam a este modelo)

model MF_ST
    // Parametros dos Degraus
    parameter Real h3sp_bias = 12.  , h3sp_degrau =  1.  ;
    parameter Real h4sp_bias = 12.  , h4sp_degrau = -1.  ;
    parameter Real h2sp_bias =  9.  , h2sp_degrau =  1.  ;   // v12: era 12 (T1/T2 a 12 cm pedem ~39 cm³/s por ramo, o máximo da B-01)
    parameter Real h1sp_bias =  9.  , h1sp_degrau = -1.  ;   // v12: era 12
    parameter Real TimeYset  = 300.  , TimeUcarga = 1000.; // TimeUcarga: usado só pelos critérios ITAE (comentados abaixo). v12: TimeYset era 10 s (degrau durante o enchimento); em 300 s os tanques já estão perto do SP
    parameter Real r1=0.15, r2 =0.15;  // v12: eram 0,6 e 0,4 (r1 + r2 = 1 deixa h1 e h2 dependentes; T5/T6 >= 5 cm pede r <= ~0,3)
    // v11: removido DIRETO (só era usado nas equações de perturbação de carga comentadas, também removidas)

    parameter Real pareamento_cruzado = 1; // v12: 1: C1 -> FV-02 e C2 -> FV-01 (indicado pelo RGA com r1, r2 < 0,5); 0: C1 -> FV-01 e C2 -> FV-02, como até a v11


    // Parametros das Bombas
      parameter Real rot_bomba_inicial_bias = 1.0;  // v12: era 0,4 (a 40% a B-01 não enche T1/T2 nem com FV-01/FV-02 abertas)
      parameter Real rot_bomba3_bias = 1.0;
      parameter Real rot_bomba4_bias = 1.0;
      parameter Real divisao_bomba_bias = 0.5;      // 50% do tempo para cada ramo

    // v12: sintonia SIMC (PI) a partir do modelo linearizado (linearizacao/analise.py); até a v11: Kp = ±0,05, Ti = 120, Td = 1
    PID_ISA C1(Yset_bias = h1sp_bias, U_bias = 0.5, Kp =  0.062, Ti = 33, Td = 0),  // nível de T1 (FV-02 no pareamento cruzado): abrir a válvula enche o tanque (Kc > 0)
            C2(Yset_bias = h2sp_bias, U_bias = 0.5, Kp =  0.062, Ti = 33, Td = 0),  // nível de T2 (FV-01 no pareamento cruzado): abrir a válvula enche o tanque (Kc > 0)
            C3(Yset_bias = h3sp_bias, U_bias = 0.5, Kp = -0.051, Ti = 40, Td = 0),  // tanque intermediario - nível sobe, tem que abrir a válvula. Pela definição do erro (Kc < 0)
            C4(Yset_bias = h4sp_bias, U_bias = 0.5, Kp = -0.051, Ti = 40, Td = 0);  // tanque intermediario - nível sobe, tem que abrir a válvula.  Pela definição do erro (Kc < 0)

    SixTanks P(r1=r1,r2=r2);

    // Criterios de Desempenho -- ignorar por enquanto
    //Real ITAE_Servo1 (start=0), ITAE_Regulatorio1 (start=0) ;
    //Real ITAE_Servo2 (start=0), ITAE_Regulatorio2 (start=0) ;

    Real e1,e2, e3, e4; // erros dos controladores
    Real h1,h2,h3,h4,h5,h6;


equation

    // -------------------------------------------------------
    // Setpoints e ligação dos blocos
    // -------------------------------------------------------

    // quatro malhas SISO (v11: o comentário da v10 dizia "SÓ LADO IMPAR"):
    //   C1 (h1) -> FV-02 (entrada, lado par)     C2 (h2) -> FV-01 (entrada, lado ímpar)   v12: pareamento cruzado
    //   C3 (h3) -> FV-03 (saída de T3)           C4 (h4) -> FV-04 (saída de T4)
    // v12: com r1, r2 < 0,5 quase toda a vazão de FV-02 chega a T1 pelo caminho T5 -> T3 -> T1
    //      (e a de FV-01 a T2 por T6 -> T4 -> T2). Com pareamento_cruzado = 0 volta o pareamento da v11.

    h1 =P.h1out; h2 =P.h2out; h3 =P.h3out;
    h4 =P.h4out; h5 =P.h5out; h6 =P.h6out;

    // tanques inferiores

    // --- controle do tanque 1 ---
    C1.SP    = h1sp_bias     + (if time < TimeYset   then 0 else h1sp_degrau);
    C1.PV    = h1;
    e1 = C1.SP - C1.PV;

    P.x1 = if pareamento_cruzado > 0.5 then C2.MV else C1.MV;  // v12: FV-01 (entrada da esquerda) com C2 no pareamento cruzado, C1 no direto


    // --- controle do tanque 2 ---
    C2.SP    = h2sp_bias     + (if time < TimeYset   then 0 else h2sp_degrau);
    C2.PV    = h2;
    e2 = C2.SP - C2.PV;

    P.x2 = if pareamento_cruzado > 0.5 then C1.MV else C2.MV;  // v12: FV-02 (entrada da direita) com C1 no pareamento cruzado, C2 no direto


    // controle tanques intermed

    // --- controle do tanque 3 ---
    C3.SP    = h3sp_bias     + (if time < TimeYset   then 0 else h3sp_degrau);

    C3.PV    = h3;
    e3 = C3.SP - C3.PV;

    P.x3 = C3.MV;  // PID abre e fecha a FV-03 (saída de T3)


    // controle do tanque 4
    C4.SP    = h4sp_bias     + (if time < TimeYset   then 0 else h4sp_degrau);
    C4.PV    = h4;
    e4 = C4.SP - C4.PV;

    P.x4 = C4.MV;  // PID abre e fecha a FV-04 (saída de T4)


    // entradas fixas
    P.rot_bomba_inicial = rot_bomba_inicial_bias;
    P.rot_bomba3 = rot_bomba3_bias;
    P.rot_bomba4 = rot_bomba4_bias;
    P.divisao_bomba_inicial = divisao_bomba_bias;
    // por enquanto travar a potencia das bombas intermed


    //der(ITAE_Servo1)= (if (time>=TimeYset and time<TimeUcarga) then abs(e1)*(time-TimeYset) else 0);
    //der(ITAE_Servo2)= (if (time>=TimeYset and time<TimeUcarga) then abs(e2)*(time-TimeYset) else 0);
    //der(ITAE_Regulatorio1)= (if (time>=TimeUcarga) then abs(e1)*(time-TimeUcarga) else 0);
    //der(ITAE_Regulatorio2)= (if (time>=TimeUcarga) then abs(e2)*(time-TimeUcarga) else 0);

end MF_ST;


// =====================================================================
// v12: planta em malha aberta, parada no ponto de operação, para linearize()
// ---------------------------------------------------------------------
// É a mesma SixTanks, sem controladores. A inicialização procura as aberturas
// x1_op ... x4_op e os níveis de T5/T6 que deixam todos os tanques em
// equilíbrio com h1 ... h4 nos setpoints. As entradas u1 ... u4 são desvios
// das aberturas em relação a esse ponto. Com u = 0 a planta fica parada, então
// linearize() dá A, B, C, D exatamente no ponto de operação.
// Os valores padrão são os do MF_ST (SP antes do degrau). Uso: linearizacao/analise.py.
// =====================================================================
model SixTanks_MA
    parameter Real r1 = 0.15, r2 = 0.15;
    parameter Real rot_bomba_inicial = 1.0, rot_bomba3 = 1.0, rot_bomba4 = 1.0;
    parameter Real divisao_bomba_inicial = 0.5;
    parameter Real h1_op = 9, h2_op = 9, h3_op = 12, h4_op = 12;  // níveis do ponto de operação (cm)

    // calculados na inicialização (fixed = false)
    parameter Real x1_op(fixed=false, start=0.5), x2_op(fixed=false, start=0.5);
    parameter Real x3_op(fixed=false, start=0.2), x4_op(fixed=false, start=0.2);

    input  Real u1(start=0), u2(start=0), u3(start=0), u4(start=0);  // desvios das aberturas de FV-01 ... FV-04
    output Real y1, y2, y3, y4, y5, y6;                              // níveis h1 ... h6 (cm)

    SixTanks P(r1=r1, r2=r2,
               h1(start=h1_op), h2(start=h2_op), h3(start=h3_op), h4(start=h4_op),
               h5(start=5, fixed=false), h6(start=5, fixed=false));

initial equation
    der(P.h1) = 0; der(P.h2) = 0; der(P.h3) = 0;
    der(P.h4) = 0; der(P.h5) = 0; der(P.h6) = 0;

equation
    P.x1 = x1_op + u1;  P.x2 = x2_op + u2;
    P.x3 = x3_op + u3;  P.x4 = x4_op + u4;
    P.rot_bomba_inicial = rot_bomba_inicial;
    P.rot_bomba3 = rot_bomba3;
    P.rot_bomba4 = rot_bomba4;
    P.divisao_bomba_inicial = divisao_bomba_inicial;

    y1 = P.h1out; y2 = P.h2out; y3 = P.h3out;
    y4 = P.h4out; y5 = P.h5out; y6 = P.h6out;
end SixTanks_MA;
