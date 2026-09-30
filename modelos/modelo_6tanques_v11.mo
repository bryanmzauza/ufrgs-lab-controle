// =====================================================================
// Modelo seis tanques -- versão 11
// ---------------------------------------------------------------------
// Mudanças em relação à v10.
// As linhas alteradas estão marcadas com "v11:".
//
//  1. Tanque vazio não esvazia mais. Na v10 a descarga era
//     CD*max(h, 1e-4)^p, que nunca chega a zero: um tanque vazio ainda
//     perdia ~0,035 cm³/s e h1/h2 ficavam negativos logo depois da partida
//     (até -0,01 cm). Agora a descarga é a função vazao_descarga, que vale
//     zero com h <= 0. Acima de 1e-3 cm ela é idêntica à da v10.
//     As descargas viraram variáveis (q1 ... q6), então aparecem no resultado.
//
//  2. B-03 e B-04 perdem sucção quando T3/T4 esvaziam. Na v10 a vazão da
//     bomba não dependia do nível, e a bomba podia tirar água de um tanque
//     vazio (nível negativo). Agora as vazões de B-03/B-04 são multiplicadas
//     por fator_succao, que vai de 0 (tanque vazio) a 1 (nível >= h_succao_min).
//
//  3. Bomba parada não bombeia. Na v10, sqrt(max(1e-6, dP)) valia 0,001 com
//     dP = 0, e sobrava vazão com rotação 0 (~0,01 cm³/s). Agora é raiz_reg(dP),
//     que vai a zero com dP = 0 e é idêntica a sqrt(dP) para dP >= 1e-6.
//
//  4. Limpeza:
//     - removidos CD3, CD4, pot_descarga3 e pot_descarga4 (não entravam em
//       nenhuma equação: T3 e T4 são esvaziados pelas bombas);
//     - removido DIRETO e as equações comentadas de perturbação de carga que
//       usavam P.F1/P.F2 como entrada (estão na v10);
//     - F_max_bomba_saida3/4 viraram parameter (eram Real com valor fixo);
//     - delta_P_max_esq e delta_P_max_dir tinham a mesma equação: viraram um só
//       delta_P_max_inicial;
//     - comentários desatualizados corrigidos ("SÓ LADO IMPAR", dicionário S1).
//
//  Não mudou: o bloco PID_ISA inteiro, equações dos balanços, curvas das bombas,
//  válvulas, valores dos parâmetros, sintonia, setpoints e ligações do MF_ST. Os estados estacionários
//  são os mesmos da v10; só os primeiros instantes da partida mudam.
//
//  Observações que continuam valendo (não são correções):
//  - PID_ISA: Id0 = Yset_bias supõe o nível já no setpoint. Com os tanques
//    partindo vazios, a derivada dá um pico Ud(0) = Kp*N*SP na partida
//    (dura ~Td/N = 0,1 s). Mantido de propósito, igual à v10.
//  - CD = 0.20*sqrt(2*981) vem de Torricelli (expoente 0,5), mas é usado com
//    expoente 0,6: na prática é um coeficiente empírico, a ajustar com a bancada.
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
    parameter Real h2sp_bias = 12.  , h2sp_degrau =  1.  ;
    parameter Real h1sp_bias = 12.  , h1sp_degrau = -1.  ;
    parameter Real TimeYset  = 10.   , TimeUcarga = 1000.; // TimeUcarga: usado só pelos critérios ITAE (comentados abaixo)
    parameter Real r1=0.6, r2 =0.4;
    // v11: removido DIRETO (só era usado nas equações de perturbação de carga comentadas, também removidas)


    // Parametros das Bombas
      parameter Real rot_bomba_inicial_bias = 0.4;  // inicio, bomba 40%
      parameter Real rot_bomba3_bias = 1.0;
      parameter Real rot_bomba4_bias = 1.0;
      parameter Real divisao_bomba_bias = 0.5;      // 50% do tempo para cada ramo

    PID_ISA C1(Yset_bias = h1sp_bias, U_bias = 0.5, Kp = 0.05, Ti = 120),  // tanque de entrada - nível sobe, tem que fechar a válvula.  Pela definição do erro (Kc > 0)
            C2(Yset_bias = h2sp_bias, U_bias = 0.5, Kp = 0.05, Ti = 120),  // tanque de entrada - nível sobe, tem que fechar a válvula. Pela definição do erro (Kc > 0)
            C3(Yset_bias = h3sp_bias, U_bias = 0.5, Kp = -0.05, Ti = 120), // tanque intermediario - nível sobe, tem que abrir a válvula. Pela definição do erro (Kc < 0)
            C4(Yset_bias = h4sp_bias, U_bias = 0.5, Kp = -0.05, Ti = 120); // tanque intermediario - nível sobe, tem que abrir a válvula.  Pela definição do erro (Kc < 0)

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
    //   C1 -> FV-01 (entrada, lado ímpar)   C2 -> FV-02 (entrada, lado par)
    //   C3 -> FV-03 (saída de T3)           C4 -> FV-04 (saída de T4)

    h1 =P.h1out; h2 =P.h2out; h3 =P.h3out;
    h4 =P.h4out; h5 =P.h5out; h6 =P.h6out;

    // tanques inferiores

    // --- controle do tanque 1 ---
    C1.SP    = h1sp_bias     + (if time < TimeYset   then 0 else h1sp_degrau);
    C1.PV    = h1;
    e1 = C1.SP - C1.PV;

    P.x1 = C1.MV;  // PID abre e fecha a FV-01 (entrada da esquerda)


    // --- controle do tanque 2 ---
    C2.SP    = h2sp_bias     + (if time < TimeYset   then 0 else h2sp_degrau);
    C2.PV    = h2;
    e2 = C2.SP - C2.PV;

    P.x2 = C2.MV;  // PID abre e fecha a FV-02 (entrada da direita)


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
