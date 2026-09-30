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
    parameter Real CD1=0.20*sqrt(2*981), CD2=0.20*sqrt(2*981),
                   CD3=0.15*sqrt(2*981), CD4=0.15*sqrt(2*981), // coeficientes maximos, vai variar conforme a abertura da válvula
				   CD5=0.20*sqrt(2*981), CD6=0.20*sqrt(2*981);

    parameter Real A1 = 63, A2 = 63, A3 = 63, A4 =63, A5=63, A6=63;  // cm² 

    parameter Real F_max_ramo = 83.0; // vazao máxima que a bomba consegue fornecer para cada um dos lados (entrada de água no sistema) (5L/min -> cm³/s)
              Real F_max_bomba_saida3 = 83.0; // vazao máxima que a bomba consegue puxar (saída de água tanque 3)   (5L/min -> cm³/s)
              Real F_max_bomba_saida4 = 83.0; // vazao máxima que a bomba consegue puxar (saída de água tanque 4)   (5L/min -> cm³/s)


    parameter Real pot_descarga1 = 0.6, pot_descarga2 = 0.6,    // potencia de descarga dos tanques. não será 0.5 pq não e totalmente turbulento
                    pot_descarga3 = 0.6, pot_descarga4 = 0.6,   // tb não será 1 pq não é totalmente laminar. 0.5 <= potencia <= 1.0 
                    pot_descarga5 = 0.6, pot_descarga6 = 0.6;


    // Variáveis do sistema
    Real F_bomba_inicial(start=10.0), F1(start=5.0), F2(start=5.0);
    Real F1_completo(start=10.0), F2_completo(start=10.0);
    Real F3(start=5.0), F4(start=5.0);
    Real F1_direto, F1_indireto, F2_direto, F2_indireto;
    Real F3_valvula_controle(start=4.0), F3_refluxo(start=1.0);
    Real F4_valvula_controle(start=4.0), F4_refluxo(start=1.0);

    Real funcao_obturador_valvula_1, funcao_obturador_valvula_2;
    Real funcao_obturador_valvula_3, funcao_obturador_valvula_4;

    Real delta_P_max_esq(start=50.0);
    Real delta_P_max_dir(start=50.0);
    Real delta_P_bomba_lado_esq(start=40.0), delta_P_bomba_lado_dir(start=40.0);
    Real delta_P_max_3(start=50.0);
    Real delta_P_max_4(start=50.0);
    Real delta_P_bomba_3(start=40.0), delta_P_bomba_4(start=40.0);
    Real F_max_atual, F_max_atual_bomba3, F_max_atual_bomba4;


    Real h1(start=0, fixed=true), h2(start=0, fixed=true),     // nível inicial em cm (0 = tanque vazio)
         h3(start=0, fixed=true), h4(start=0, fixed=true),     // nível inicial em cm (0 = tanque vazio)
         h5(start=0, fixed=true), h6(start=0, fixed=true);     // nível inicial em cm (0 = tanque vazio)

    
    // criando variáveis para atuar -- abertura das válvulas
    input Real x1; // Abertura da válvula após a bomba (Perna Esquerda)
    input Real x2; // Abertura da válvula após a bomba (Perna Direita)
    input Real x3; // Abertura da válvula de saída do T3
    input Real x4; // Abertura da válvula de saída do T4

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

    // BOMBA 01 - PUXA ÁGUA DO RESERVATÓRIO

    // ideia será calcular qual a vazão em cada ramo se fosse só por aquele ramo e ponderar por quanto tempo ele vai para cada ramo
    // variável "completo" é caso fosse tudo para aquele lado


    // lado esquerdo
    delta_P_max_esq = K_bomba_inicial *rot_bomba_inicial^2;  // Delta P max é função da rotação da bomba 
    delta_P_bomba_lado_esq = delta_P_max_esq * sqrt(max(1e-6, 1.0 - (F1_completo / max(1e-5, F_max_atual))^2));                                              // eq da bomba
    F1_completo = funcao_obturador_valvula_1 * CV_valvula1 * sqrt(max(1e-6, delta_P_bomba_lado_esq));


    // lado direito

    delta_P_max_dir = K_bomba_inicial * rot_bomba_inicial^2;  // Delta P max é função da rotação da bomba 
    delta_P_bomba_lado_dir = delta_P_max_dir * sqrt(max(1e-6, 1.0 - (F2_completo / max(1e-5, F_max_atual))^2));                                              // eq da bomba
    F2_completo = funcao_obturador_valvula_2 * CV_valvula2 * sqrt(max(1e-6, delta_P_bomba_lado_dir));

    F1 = F1_completo * divisao_bomba_inicial;       // ponderação pelo tempo que joga água para a esq
    F2 = F2_completo * (1 - divisao_bomba_inicial); // ponderação pelo tempo que joga água para a dir

    F_bomba_inicial = F1 + F2; // conservação de massa


    // --- Divisão Válvulas de 3 Vias --- r é o quanto vai para o tanque direto
    F1_direto = r1 * F1;          // quanto vai de F1 para o tanque 3
    F1_indireto = (1 - r1) * F1;  // quanto vai de F1 para o tanque 6

    F2_direto = r2 * F2;          // quanto vai de F2 para o tanque 4
    F2_indireto = (1 - r2) * F2;  // quanto vai de F2 para o tanque 5


    
	A5 * der ( h5 ) = F2_indireto - CD5*(max(h5, 1e-4)^pot_descarga5);  // descarga por gravidade
    A6 * der ( h6 ) = F1_indireto - CD6*(max(h6, 1e-4)^pot_descarga6);  // descarga por gravidade



    // tanques intermediarios - tem a valvula intermediaria e bomba

    // BOMBAS 03 E 04 -- INTERMEDIARIAS

    F_max_atual_bomba3 = rot_bomba3 * F_max_bomba_saida3; // considerando, inicialmente, que a vazão máxima varia linearmente com a rotação da bomba 
    F_max_atual_bomba4 = rot_bomba4 * F_max_bomba_saida4; // considerando, inicialmente, que a vazão máxima varia linearmente com a rotação da bomba 

    funcao_obturador_valvula_3 = max(0.0, x3); // obturador linear por enquanto
    funcao_obturador_valvula_4 = max(0.0, x4); // obturador linear por enquanto
    

    delta_P_max_3 = K_bomba_3 * rot_bomba3^2;    // Delta P max é função da rotação da bomba 
    delta_P_bomba_3 = delta_P_max_3 * sqrt(max(1e-6, 1.0 - (F3 / max(1e-5, F_max_atual_bomba3))^2));             // curva da bomba

    delta_P_max_4 = K_bomba_4 * rot_bomba4^2;    // Delta P max é função da rotação da bomba 
    delta_P_bomba_4 = delta_P_max_4 * sqrt(max(1e-6, 1.0 - (F4 / max(1e-5, F_max_atual_bomba4))^2));             // curva da bomba

    // Equacoes de saída

    F3 = F3_refluxo + F3_valvula_controle;
    F3_valvula_controle = funcao_obturador_valvula_3 * CV_valvula3 * sqrt(max(1e-6, delta_P_bomba_3));
    F3_refluxo = CV_valvula3_refluxo * abertura_valv_refluxo_3 * sqrt(max(1e-6, delta_P_bomba_3));

    
    F4 = F4_refluxo + F4_valvula_controle;
    F4_valvula_controle = funcao_obturador_valvula_4 * CV_valvula4 * sqrt(max(1e-6, delta_P_bomba_4));
    F4_refluxo = CV_valvula4_refluxo * abertura_valv_refluxo_4 * sqrt(max(1e-6, delta_P_bomba_4));


    // BALANÇO TANQUES INTERMEDIARIOS E SAÍDA

    A4 * der ( h4 ) = CD6*(max(h6, 1e-4)^pot_descarga6) + F2_direto - F4_valvula_controle; // entra o que sai do tanque 6, o desvio do refluxo de F2 (entrada) e o refluxo e sai o que a bomba puxa
    A3 * der ( h3 ) = CD5*(max(h5, 1e-4)^pot_descarga5) + F1_direto - F3_valvula_controle; // entra o que sai do tanque 5, o desvio do refluxo de F1 (entrada) e o refluxo e sai o que a bomba puxa. Notar que o refluxo sai, mas entra de volta. Logo, ao somar F3_refluxo e subtrair F_saida3 = F_refluxo + F3_sai, tem-se apenas F3_sai...

	A2 * der ( h2 ) = F4_valvula_controle - CD2*(max(h2, 1e-4)^pot_descarga2);
    A1 * der ( h1 ) = F3_valvula_controle - CD1*(max(h1, 1e-4)^pot_descarga1);

    h1out = h1; h2out = h2; h3out = h3;
    h4out = h4; h5out = h5; h6out = h6;
    
end SixTanks;
//S1 = {'TimeYset': 5.,'h1sp_degrau':1.,'h2sp_degrau':-1,'TimeUcarga':100.,'F1_degrau':1.,'F2_degrau': -1, 'h1sp_bias':12, 'h2sp_bias':12,'r1':0.6,'r2':0.4,'Direct':true}

model MF_ST
    // Parametros dos Degraus
    parameter Real h3sp_bias = 12.  , h3sp_degrau =  1.  ;
    parameter Real h4sp_bias = 12.  , h4sp_degrau = -1.  ;
    parameter Real h2sp_bias = 12.  , h2sp_degrau =  1.  ;
    parameter Real h1sp_bias = 12.  , h1sp_degrau = -1.  ;
    parameter Real TimeYset  = 10.   , TimeUcarga = 1000.;
    //parameter Real F1_bias   = 10.    ,  F1_degrau = 1.;
    //parameter Real F2_bias   = 10.    ,  F2_degrau = 1.;
    parameter Real r1=0.6, r2 =0.4;
    parameter Real DIRETO = 1.;   


    // Parametros das Bombas
      parameter Real rot_bomba_inicial_bias = 0.4;  // incio, bomba 40%
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
    // Perturbacoes (Setpoint e Carga) e conectando os blocos
    // -------------------------------------------------------
    // P.r1 = r1;
    // P.r2 = r2;


    // fazendo SISO - SÓ LADO IMPAR
       
    h1 =P.h1out; h2 =P.h2out; h3 =P.h3out;
    h4 =P.h4out; h5 =P.h5out; h6 =P.h6out;

    // tanques inferiores
    
    // --- controle do tanque 1 ---
    C1.SP    = h1sp_bias     + (if time < TimeYset   then 0 else h1sp_degrau);
    C1.PV    = h1;
    e1 = C1.SP - C1.PV;
            
    P.x1 = C1.MV;  // PID abre e fecha a válvula entrada da esquerda


    // --- controle do tanque 2 ---
    C2.SP    = h2sp_bias     + (if time < TimeYset   then 0 else h2sp_degrau);
    C2.PV    = h2;
    e2 = C2.SP - C2.PV;
            
    P.x2 = C2.MV;  // PID abre e fecha a válvula entrada da dir
    

    // controle tanques intermed

    // --- controle do tanque 3 ---
    C3.SP    = h3sp_bias     + (if time < TimeYset   then 0 else h3sp_degrau);

    C3.PV    = h3;
    e3 = C3.SP - C3.PV;
            
    P.x3 = C3.MV;  // PID abre e fecha a válvula esquerda


    // controle do tanque 4
    C4.SP    = h4sp_bias     + (if time < TimeYset   then 0 else h4sp_degrau);
    C4.PV    = h4;
    e4 = C4.SP - C4.PV;
            
    P.x4 = C4.MV;  // PID abre e fecha a válvula dir
    



    // P.F1  + (if time < TimeUcarga then 0 else F1_degrau) = (if DIRETO > 0 then C1.MV else C2.MV);        
    // P.F2  + (if time < TimeUcarga then 0 else F2_degrau) = (if DIRETO > 0 then C2.MV else C1.MV);
   
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



