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
    
    // Process parameters - parâmetros das descargas dos tanques
    parameter Real CD1=0.20*sqrt(2*981), CD2=0.20*sqrt(2*981),
                   CD3=0.15*sqrt(2*981), CD4=0.15*sqrt(2*981), // coeficientes maximos, vai variar conforme a abertura da válvula
				   CD5=0.20*sqrt(2*981), CD6=0.20*sqrt(2*981);

    parameter Real A1 = 63, A2 = 63, A3 = 63, A4 =63, A5=63, A6=63;  // cm² 

    parameter Real F_max_ramo = 83.0; // vazao máxima que a bomba consegue fornecer para cada um dos lados (entrada de água no sistema)
              Real F_max_bomba_saida3 = 83.0; // vazao máxima que a bomba consegue puxar (saída de água tanque 3)
              Real F_max_bomba_saida4 = 83.0; // vazao máxima que a bomba consegue puxar (saída de água tanque 4)

    parameter Real pot_descarga1 = 0.5, pot_descarga2 = 0.5,    // potencia de descarga dos tanques. não será 0.5 pq não e totalmente turbulento
                    pot_descarga3 = 0.5, pot_descarga4 = 0.5,   // tb não será 1 pq não é totalmente laminar. 0.5 <= potencia <= 1.0 
                    pot_descarga5 = 0.5, pot_descarga6 = 0.5;



// modelar os tanques intermediarios -- $$P_{bomba} = P_{max} \cdot (\text{pot\_bomba})^2 - K_b \cdot F_{total}^2$$
    parameter Real P_max = 80.0;    // Pressão máxima da bomba a 100% de potência (cmH2O)
    parameter Real K_b = 0.02;      // Coeficiente de perda interna da bomba
    parameter Real K_ref3 = 3.0;     // Condutância do refluxo tanque3 (relacionado com a abertura da válvula que tem no refluxo)
    parameter Real K_ref4 = 3.0;     // Condutância do refluxo tanque4 (relacionado com a abertura da válvula que tem no refluxo)
    Real F_saida4, F_saida3, F3_sai, F4_sai, F3_refluxo, F4_refluxo; // variaveis relacionadas com a saída dos tanques intermed




    Real h1 (start=10), h2(start=10), h3 (start =10), h4 (start=10), h5(start=10),h6(start=10); // inicializando alturas
    Real F1, F2; // incializando vazoes de cada lado do sistema

    
    input Real v1; // Abertura da válvula após a bomba (Perna Esquerda)
    input Real v2; // Abertura da válvula após a bomba (Perna Direita)
    input Real v3; // Abertura da válvula de saída do T3
    input Real v4; // Abertura da válvula de saída do T4

    input   Real r1, r2;    // divisão das correntes de entrada

    input Real pot_bomba_inicial;       // 1 = ligado, 0 = desligado, 0.5 intermediario -- input para poder usar como var de controle
    input Real pot_bomba3, pot_bomba4;  // 1 = ligado, 0 = desligado, 0.5 intermediario -- input para poder usar como var de controle

    output  Real h1out, h2out, h3out, h4out, h5out, h6out;  // saidas dos níveis

equation
    F1 = F_max_ramo * pot_bomba_inicial * v1;   // quanto a parte esq consegue puxar do reservatório principal
    F2 = F_max_ramo * pot_bomba_inicial * v2;   // quanto a parte dir consegue puxar do reservatório principal


	A5 * der ( h5 ) = F2*r2 - CD5*(max(h5,0)^pot_descarga5);
    A6 * der ( h6 ) = F1*r1 - CD6*(max(h6,0)^pot_descarga6);

    // tanques intermediarios - tem a valvula intermediaria e bomba
    F_saida4 = if h4 > 0 then (F_max_bomba_saida4 * pot_bomba4) else 0; // o que sai é o que a bomba consegue puxar * potencia que escolhi (0 a 1)  -- daria para colocar curva da bomba...
    F_saida3 = if h3 > 0 then (F_max_bomba_saida3 * pot_bomba3) else 0; // o que sai é o que a bomba consegue puxar * potencia que escolhi (0 a 1)  -- daria para colocar curva da bomba...
    
    F4_sai = F_saida4 * (v4 *CD4)/(v4*CD4 + K_ref4); // balanço de resistencia na linha -> o que segue adiante no processo é a ressitência do que segue adiante / resistencia total
    F4_refluxo = F_saida4 * (K_ref4)/(v4*CD4 + K_ref4);


    F3_sai = F_saida3 * (v3 *CD3)/(v3*CD3 + K_ref3); // balanço de resistencia na linha -> o que segue adiante no processo é a ressitência do que segue adiante / resistencia total
    F3_refluxo = F_saida3 * (K_ref3)/(v3*CD3 + K_ref3);


    A4 * der ( h4 ) = CD6*(max(h6,0)^pot_descarga6) + (1-r2)*F2 + F4_refluxo - F_saida4; // entra o que sai do tanque 6, o desvio do refluxo de F2 (entrada) e o refluxo e sai o que a bomba puxa
    A3 * der ( h3 ) = CD5*(max(h5,0)^pot_descarga5) + (1-r1)*F1 + F3_refluxo - F_saida3; // entra o que sai do tanque 5, o desvio do refluxo de F1 (entrada) e o refluxo e sai o que a bomba puxa. Notar que o refluxo sai, mas entra de volta. Logo, ao somar F3_refluxo e subtrair F_saida3 = F_refluxo + F3_sai, tem-se apenas F3_sai...

	A2 * der ( h2 ) = F4_sai - CD2*(max(h2,0)^pot_descarga2);
    A1 * der ( h1 ) = F3_sai - CD1*(max(h1,0)^pot_descarga1);

    h1out = h1;
    h2out = h2;
	h3out = h3;
    h4out = h4;
	h5out = h5;
    h6out = h6;
    
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

    // parametros da bomba inicial e valvulas inciais
    parameter Real pot_bomba_inicial_bias = 0.4;    // incio, bomba 40%
    parameter Real v1_bias = 1.0;                   // valv aberta
    parameter Real v2_bias = 1.0;                   // valv aberta

    // pot bombas intermed
    parameter Real pot_bomba3_bias = 1.0;   // comeca na potencia maxima
    parameter Real pot_bomba4_bias = 1.0;   // comeca na potencia maxima
    
    PID_ISA C1 (Yset_bias = h3sp_bias, U_bias = 0.5), // bias 0.5 para ser valvula meio aberta - controlador nivel 3
            C2 (Yset_bias = h4sp_bias, U_bias = 0.5), // bias 0.5 para ser valvula meio aberta - controlador nivel 4
            C3 (Yset_bias = h1sp_bias, U_bias = 0.5), // bias 0.5 para ser valvula meio aberta - controlador nivel 1
            C4 (Yset_bias = h2sp_bias, U_bias = 0.5); // bias 0.5 para ser valvula meio aberta - controlador nivel 2
        
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
       
    h1 =P.h1out;
    h2 =P.h2out;
    h3 =P.h3out;
    h4 =P.h4out;
    h5 =P.h5out;
    h6 =P.h6out;

    
    // controle do tanque 3
    C1.SP    = h3sp_bias     + (if time < TimeYset   then 0 else h3sp_degrau);
    C1.PV    = h3;
    e1 = C1.SP - C1.PV;
            
    P.v3 = C1.MV;  // PID abre e fecha a válvula da esquerda


    // controle do tanque 4
    C2.SP    = h4sp_bias     + (if time < TimeYset   then 0 else h4sp_degrau);
    C2.PV    = h4;
    e2 = C2.SP - C2.PV;
            
    P.v4 = C2.MV;  // PID abre e fecha a válvula da dir
    

    // controle tanques inferiores

    // controle do tanque 1
    C3.SP    = h1sp_bias     + (if time < TimeYset   then 0 else h1sp_degrau);

    C3.PV    = h1;
    e3 = C3.SP - C3.PV;
            
    P.v1 = C3.MV;  // PID abre e fecha a válvula de entrada esquerda


    // controle do tanque 2
    C4.SP    = h2sp_bias     + (if time < TimeYset   then 0 else h2sp_degrau);
    C4.PV    = h2;
    e4 = C4.SP - C4.PV;
            
    P.v2 = C4.MV;  // PID abre e fecha a válvula de entrada dir
    



    // P.F1  + (if time < TimeUcarga then 0 else F1_degrau) = (if DIRETO > 0 then C1.MV else C2.MV);        
    // P.F2  + (if time < TimeUcarga then 0 else F2_degrau) = (if DIRETO > 0 then C2.MV else C1.MV);
   
    // travar outras variaveis
    P.pot_bomba_inicial = pot_bomba_inicial_bias;  
    // P.v2 = v2_bias;         // Válvula dir fixa
    // P.v1 = v1_bias;         // Válvula esq fixa
    P.pot_bomba4 = pot_bomba4_bias;  // potencia alta para ter rangeabilidade...
    P.pot_bomba3 = pot_bomba3_bias;  // potencia alta para ter rangeabilidade...
// por enquanto travar a potencia das bombas intermed



    //der(ITAE_Servo1)= (if (time>=TimeYset and time<TimeUcarga) then abs(e1)*(time-TimeYset) else 0); 
    //der(ITAE_Servo2)= (if (time>=TimeYset and time<TimeUcarga) then abs(e2)*(time-TimeYset) else 0); 
    //der(ITAE_Regulatorio1)= (if (time>=TimeUcarga) then abs(e1)*(time-TimeUcarga) else 0);
    //der(ITAE_Regulatorio2)= (if (time>=TimeUcarga) then abs(e2)*(time-TimeUcarga) else 0);

end MF_ST;



