% il paziente 2 migliora da posizione 9
% il paziente 3 cala dalla 4 e migliora dalla 
% il paziente 4 talvolta migliora dalla posizione 5

clear all
close all
clc

global V1 V2 V3 t_dop y i ke3 k31 alpha beta gamma
global c1 c3_delay tfreq Tapping_f


alpha = 0.75;  %(0.2*(Ugo_trigger-0.8)+0.5)/(0.7*(Ugo_trigger-0.8));

%guadagno da DA a No-Go (inibizione)
beta = -1;
%guadagno da DA a interneurone colinergico (inibizione)
gamma = -0.5; 
Ns = 4;
S = zeros(Ns,1);
S(1) = 1;

paziente = 4;

switch paziente
       
    case 1   % paziente stabile
    load('parametri_pazienteG1S5_max.mat')
    %dopamina tonica
   %Dop_tonic_PD = 0.43;
   z=5; %posizione tra le 10 prove random in cui si ha la funzione costo migliore
   
    case 2
    load('parametri_pazienteG1S9_max.mat')
    %dopamina tonica
   %Dop_tonic_PD = 0.43;
   z=8; %posizione tra le 10 prove random in cui si ha la funzione costo migliore
   
    case 3
    load('parametri_pazienteG1S10_max.mat')
    %dopamina tonica
   %Dop_tonic_PD = 0.43;
   z=9; %posizione tra le 10 prove random in cui si ha la funzione costo migliore
    
    case 4   % paziente instabile
    load('parametri_pazienteG2S5_max.mat')
    %per il paziente G2s5 z= 10
    z=10; %posizione tra le 10 prove random in cui si ha la funzione costo migliore
    
    case 5
    load('parametri_pazienteG2S12_max.mat')
    %per il paziente G2s12 z= 7
    z=7; %posizione tra le 10 prove random in cui si ha la funzione costo migliore
     
end
    
    Delay=p1_totale(z);
    ke3= p2_totale(z);
    Dop_max = p3_totale(z);
    Dop_50 = p4_totale(z);
    N = p5_totale(z);
    Dop_tonic = p6_totale(z);
    


STN_ON = 1;
T_ON = 1;

ft_tot = [];
ft_tot_add = [];
DA_tot = [];

dt = 0.1;
Calculate_levodopa
tau=30;   % da togliere??
%Delay = p1_totale;
Delay_indice = round(Delay/dt);


c3_delay = [zeros(1,Delay_indice) c3(1:L-Delay_indice)']';


%inizializzazione
Nc = 4;
load W_tot_new_W0e5_D1e0
    Wgc = squeeze(Wgc_epocs(:,:,100));
    Wgs = squeeze(Wgs_epocs(:,:,100));
    Wnc = squeeze(Wnc_epocs(:,:,100));
    Wns = squeeze(Wns_epocs(:,:,100));
Ke = 7;
    
%     Dop_tonic = Dop_tonic_PD;  % devo scriverlo di nuovo perch� ho memorizzato Dop_tonic nel file appena letto


% simulo senza addestrare
for iiii=1:150:length(c3_delay)
Dop_ex = c3_delay(iiii);
% assegno dopamina tot
DA = Dop_tonic+(Dop_max*Dop_ex^N)/(Dop_50^N+Dop_ex^N);
[Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,tt,k_tap_vett,Uchi,ChI,ft] = BG_model_function_tapping_mauro_3(S,Wgc,Wgs,Wnc,Wns,Ke,STN_ON,T_ON,DA);


ft_tot = [ft_tot ft];

end

% simulo addestrando

t_addestramento= input('istant of training (min)? [it must be a multiple of 15]'); %training in minutes
i_addestramento =  150*(t_addestramento/15)+1; %index of the training;

for iiii=1:150:length(c3_delay)
Dop_ex = c3_delay(iiii);
% assegno dopamina tot
DA = Dop_tonic+(Dop_max*Dop_ex^N)/(Dop_50^N+Dop_ex^N);
[Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,tt,k_tap_vett,Uchi,ChI,ft] = BG_model_function_tapping_mauro_3(S,Wgc,Wgs,Wnc,Wns,Ke,STN_ON,T_ON,DA);   

DA_tot = [DA_tot DA];

if iiii == i_addestramento
    Synapse_training   % addestro le sinapsi a partire dai valori attuali
    Wgc = squeeze(Wgc_epocs(:,:,end));
    Wgs = squeeze(Wgs_epocs(:,:,end));
    Wnc = squeeze(Wnc_epocs(:,:,end));
    Wns = squeeze(Wns_epocs(:,:,end));    
    clear Wgc_epocs Wgc_post Wgs_epocs Wgs_post Wnc_epocs Wnc_post Wns_epocs Wns_pos;
    S = zeros(Ns,1);
    S(1) = 1;
end

ft_tot_add  = [ft_tot_add  ft];

end


t2 = t1(1:150:end);

width = 1.5;
font = 18;

close all
figure
plot(t2,ft_tot.*60,'b-o',t2,ft_tot_add.*60,'r--+','linewidth', width)
xlabel('time (min)','fontsize',font)
ylabel('tapping frequency (taps/min)','fontsize',font) 
set(gca,'fontsize',font)


