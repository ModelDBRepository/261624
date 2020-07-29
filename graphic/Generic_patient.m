
clear all
close all
clc

%% variabili globali
global V1 V2 V3 t_dop y i ke3 k31 alpha beta gamma
global S Wgc Wgs Wnc Wns Ke STN_ON T_ON 
global c3_delay tfreq Tapping_f Dop_tonic 


%% inizializzazione sinapsi
Nc = 4;
load W_tot_new_W0e5_D1e0
    Wgc = squeeze(Wgc_epocs(:,:,100));
    Wgs = squeeze(Wgs_epocs(:,:,100));
    Wnc = squeeze(Wnc_epocs(:,:,100));
    Wns = squeeze(Wns_epocs(:,:,100));
Ke = 7;


 nome_paziente = input('Name of patient file among apices? ');
 stringa1 = ['load ' nome_paziente];
 eval(stringa1)
 nome_parametri = input('Name of parameter file among apices? ');
 stringa2 = ['load ' nome_parametri];
 eval(stringa2)

%% calcolo la cinetica della levodopa
V1 = 12;
V2 = 32;
V3 = 2;
k31 = 0.02;
%ke3 = 0.03;

i = [zeros(Delay_levodopa,1);  3.33*(ones(300,1)); zeros(2200 - Delay_levodopa,1)]';

ke1 = ketot - k31;


dt = 0.1;
t1 = [0:dt:250];
L = length(t1);
c1 = zeros(L,1);   % plasma+periferico
c2 = zeros(L,1);   
c3 = zeros(L,1);
%c1(1) = D/V1;


for j = 1: L-1,
    dc1 = -k21/V1*c1(j)+k12/V1*c2(j)-ke1/V1*c1(j)-k31/V1*c1(j)+i(j)/V1;
    dc2= k21/V2*c1(j)-k12/V2*c2(j);
    c1(j+1) = c1(j) +dt*dc1;
    c2(j+1) = c2(j) + dc2*dt;
end

%index = [1 151 301 451 601 751 901 1201 1501 1801];  
index = t_dop*10+1;
kk = length(y);
index = index(1:kk);  % inserito perch? in qualche caso y ? meno lungo

%% faccio le figure della levodopa
width = 1.0;
font = 14;

figure
subplot(221)
plot(t1(index),c1(index),'bo-',t_dop,y,'r*--','linewidth', width)
%xlabel('time (min)','fontsize',font)
ylabel('\mu g/mL','fontsize',font)
title('Levodopa concentration','fontsize',font)
xlabel('time (min)','fontsize',font)
axis([0 180 0 max(y)+0.5])
set(gca,'fontsize',font)
set(gca,'XTick',[0 60 120 180])
        
M_y = mean(y,'omitnan'); %valor medio dei dati sperimentali
SS_tot = sum((y-M_y).^2,'omitnan');  % somma dei quadrati dei dati sperimentali
SS_res = sum((c1(index)-y).^2,'omitnan');  % somma dei quadrati dei residui
r2 = 1 - SS_res/SS_tot;
disp(r2)
        



%% parametri dei gangli della base
alpha = 0.75;  %(0.2*(Ugo_trigger-0.8)+0.5)/(0.7*(Ugo_trigger-0.8));

%guadagno da DA a No-Go (inibizione)
beta = -1;

%guadagno da DA a interneurone colinergico (inibizione)
gamma = -0.5;  

Ns = 4;
S = zeros(Ns,1);
S(1) = 1;
   

STN_ON = 1;
T_ON = 1;

dt = 0.1;
tau=30;


[m k1 ]= min(fcosto_totale);
        
Delay = p1_totale(k1);
ke3 = p2_totale(k1);
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

Dop_max = p3_totale(k1);
Dop_50 = p4_totale(k1);
N = p5_totale(k1);
Dop_tonic = p6_totale(k1);

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


subplot(222)
plot(tfreq,ft_tot.*60,'bo-',tfreq,Tapping_f,'r*--','linewidth', width)
%xlabel('time (min)','fontsize',font)
ylabel('taps/min','fontsize',font)
xlabel('time (min)','fontsize',font)
title('Tapping frequency','fontsize',font)
%axis([0 240 80 200])
set(gca,'fontsize',font)
set(gca,'XTick',[0 60 120 180 240])
        

M_ft = mean(Tapping_f); %valor medio dei dati sperimentali
SS_tot = sum((Tapping_f-mean(Tapping_f)).^2);  % somma dei quadrati dei dati sperimentali
SS_res = sum((60*ft_tot'-Tapping_f).^2);  % somma dei quadrati dei residui
r2 = 1 - SS_res/SS_tot;
disp(r2)
  







