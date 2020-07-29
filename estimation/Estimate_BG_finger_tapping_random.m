clear all
close all
clc

global V1 V2 V3 t_dop y i ke3 k31 alpha beta gamma
global S Wgc Wgs Wnc Wns Ke STN_ON T_ON 
global c1 c3_delay tfreq Tapping_f Dop_tonic


%% inizializzazione sinapsi
Nc = 4;
load W_tot_new_W0e5_D1e0
    Wgc = squeeze(Wgc_epocs(:,:,100));
    Wgs = squeeze(Wgs_epocs(:,:,100));
    Wnc = squeeze(Wnc_epocs(:,:,100));
    Wns = squeeze(Wns_epocs(:,:,100));
Ke = 7;

 

%% carico i dati e i parametri del paziente
nome_paziente = input('Name of the patient data file (between apices)? ');
stringa = ['load ' nome_paziente];
eval(stringa)



%% parametri basali
width = 1.5;
font = 18;

alpha = 0.75;  %(0.2*(Ugo_trigger-0.8)+0.5)/(0.7*(Ugo_trigger-0.8));

%guadagno da DA a No-Go (inibizione)
beta = -1;
%guadagno da DA a interneurone colinergico (inibizione)
gamma = -0.5; 

Ns = 4;
S = zeros(Ns,1);
S(1) = 1;

% calcolo il valore di DOp_tonic dalla curva del finger tapping
load curva_tapping_3
Tapping_f_iniziale = mean(Tapping_f(1:1));  % posso scegliere di mediare fra più valori
vett_index = find(Tapping_f_iniziale> Freq);
index = max(vett_index);
Dop_tonic0 = Valori_dopamina(index) + (Valori_dopamina(index+1) - Valori_dopamina(index) ) * ( Tapping_f_iniziale - Freq(index) ) / ( Freq(index+1) - Freq(index) );
if Dop_tonic0 < 0.41  
    Dop_tonic0 = 0.41; 
end
% valori iniziale degli altri parametri

    

STN_ON = 1;
T_ON = 1;

Np = 10;

p1_totale = zeros(Np,1);
p2_totale = zeros(Np,1);
p3_totale = zeros(Np,1);
p4_totale = zeros(Np,1);
p5_totale = zeros(Np,1);
p6_totale = zeros(Np,1);
fcosto_totale = zeros(Np,1);

dt = 0.1;

%% eseguo la stima dei parametri della concentrazione ematica di levodopa

V1 = 12;
V2 = 32;
V3 = 2;
k31 = 0.02;
ke3 = 0.03;

% parametri [Delay_levodopa k21 k12 ketot]   
% Delay_levodopa è diviso 100 % per aver un valore confrontabile con gli altri parametri
p0 = [2 1.5 1.5 3];
[p fval] = fminsearch('Cost_levodopa',p0);


%% simulo con i parametri ottenuti durante l'addestramento
Delay_levodopa = round(p(1)*100);   
k21 = p(2);
k12 = p(3);
ketot = p(4);
ke1 = ketot - k31;

i = [zeros(Delay_levodopa,1);  3.33*(ones(300,1)); zeros(2200 - Delay_levodopa,1)]';

dt = 0.1;
t1 = [0:dt:250];
L = length(t1);
c1 = zeros(L,1);   % plasma+periferico
c2 = zeros(L,1);   


for j = 1: L-1,
    dc1 = -k21/V1*c1(j)+k12/V1*c2(j)-ke1/V1*c1(j)-k31/V1*c1(j)+i(j)/V1;
    dc2= k21/V2*c1(j)-k12/V2*c2(j);
    c1(j+1) = c1(j) +dt*dc1;
    c2(j+1) = c2(j) + dc2*dt;
end

tau=30;





%% eseguo la stima dei parametri della DA
% eseguo un ciclo for variando sia N che Dop_50
% parametri [Delay Dop_max Dop_50 N Dop_tonic0]



Hfigure=0;  %gparametro della figura
p0(6) = Dop_tonic0;

for k1 = 1:Np,
p0(1) = rand(1,1)*35;
p0(2)= rand(1,1)*0.06;
p0(3) = rand(1,1)*0.3;
p0(4) = rand(1,1)*1.2;
p0(5) = rand(1,1)*12;


        
ft_tot = [];
DA_tot = [];

options = optimset('Display','iter','TolFun',1e-1,'MaxIter',300);
fun = @Cost_tapping;
[p, fval] = fminsearch(fun,p0,options);



%% simulo con i parametri ottenuti per l'addestramento 


Delay = p(1);
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
ke3 = p(2);
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

%% faccio le figure
Hfigure = Hfigure+1;
figure(Hfigure)
plot(tfreq,ft_tot.*60,'b-',tfreq,Tapping_f,'r*--','linewidth', width)
xlabel('time (min)','fontsize',font)
ylabel('tapping frequency (taps/min)','fontsize',font) 
set(gca,'fontsize',font)

p1_totale(k1) = p(1);
p2_totale(k1) = p(2);
p3_totale(k1) = p(3);
p4_totale(k1) = p(4);
p5_totale(k1) = p(5);
p6_totale(k1) = p(6);
fcosto_totale(k1) = fval;

end

save parametri_paziente Delay_levodopa k21 k12 ke3 ketot p1_totale p2_totale p3_totale p4_totale p5_totale p6_totale fcosto_totale Np


