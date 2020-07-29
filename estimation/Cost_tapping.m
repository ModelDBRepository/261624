function costo = cost_tapping(p)

global V1 V2 V3 t_dop y i ke3 k31 alpha beta gamma
global S Wgc Wgs Wnc Wns Ke STN_ON T_ON
global c1 c3_delay tfreq Tapping_f Dop_tonic
% calcola la funzione costo per la stima dei parametri della dopamina

minimo = min(p);
if minimo > 0
    Delay = p(1);
    ke3 = p(2);
    
    dt = 0.1;
    t1 = [0:dt:250];
    L = length(t1);
    c3 = zeros(L,1);
    
    for j = 1: L-1,
        dc3= k31/V3*c1(j) -ke3/V3*c3(j);
        c3(j+1) = c3(j) + dc3*dt;
    end
    
    Delay_indice = round(Delay/dt);
    
    % variabile c3_ritardata
    c3_delay = [zeros(1,Delay_indice) c3(1:L-Delay_indice)']';
    %il passo di integrazione per il modello della levodopo era 0.1. Per cui
    % ricampiono a un minuto;
    c3_delay = c3_delay(1:10:length(t1));
    
    
    Dop_max = p(3);
    Dop_50 = p(4);
    N = p(5);
    Dop_tonic = p(6);
    
    DA_tot = [];
    ft_tot = [];
    for iiii=1:1:length(tfreq)
        index = tfreq(iiii)+1;
        Dop_ex = c3_delay(index);
        % assegno dopamina tot
        DA = Dop_tonic+(Dop_max*Dop_ex^N)/(Dop_50^N+Dop_ex^N);
        [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,tt,k_tap_vett,Uchi,ChI,ft] = BG_model_function_tapping_mauro_3(S,Wgc,Wgs,Wnc,Wns,Ke,STN_ON,T_ON,DA);
        
        DA_tot = [DA_tot DA];
        ft_tot = [ft_tot ft];
        
    end
    
    diff  = (ft_tot'*60 - Tapping_f);
    M = max(abs(diff));
    
    K_max = 10;    % valore che può essere modificato
    costo = sum(diff.^2)+ K_max*M;
else
    costo = 1E8;
end
end

