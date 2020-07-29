function [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,k_tap_vett,Uchi,ChI,ft] = BG_model_function_tapping_mauro_3(S,Wgc,Wgs,Wnc,Wns,Ke,STN_ON,T_ON,Dop_tonic)

% function che simula un soggetto con dinamica pi� veloce (tesi di Alfredo
% Prete e Barbara Ferro)

global V1 V2 V3 t_dop y i ke3 k31 alpha beta gamma
global c1 c3_delay tfreq Tapping_f

tau = 10; %12;   %costante di tempo di base  15
tauL = 5*tau;   %costante di tempo dell'inibizione laterale
%tauS = tau/5;

dt = 0.1;   %passo
t = (0:dt:1000)';   %tempi [ms]
%t = (0:dt:3000)';
%t = (0:dt:180)';   %tempi [ms]
D = length(t);   %numero campioni


%% blocco di inizializzazione strutture

%C: corteccia
Nc = 4;   %neuroni della corteccia
C = zeros(Nc,D);   %uscita neuroni
Uc = zeros(Nc,D);   %ingresso a sigmoide
Ul = zeros(Nc,D);   %contributo inibizione laterale
%Go: striato, Go
Go = zeros(Nc,D);
Ugo = zeros(Nc,D);
%NoGo: striato, No-Go
NoGo = zeros(Nc,D);
Unogo = zeros(Nc,D);
%Gpe: globo pallido parte esterna
Gpe = zeros(Nc,D);
Ugpe = zeros(Nc,D);
%Gpi: globo pallido parte interna
Gpi = zeros(Nc,D);
Ugpi = zeros(Nc,D);
%T: talamo
T = zeros(Nc,D);
Ut = zeros(Nc,D);
%STN: nucleo sub-talamico
STN = zeros(1,D);
Ustn = zeros(1,D);
%E: energia (conflitto corteccia)
E = zeros(1,D);
%ChI: interneuroni colinergici
ChI = zeros(1,D);
Uchi = zeros(1,D);

%ingresso DA+ACh a Go e NoGo
IGo_DA_Ach = zeros(Nc,D);
INoGo_DA_Ach = zeros(Nc,D);


%% blocco di inizializzazione sinapsi

%pesi da stimolo a cortecciaF
Wcs = 1*ones(4,4);   %0.2 extradiag, 1.1 su diag
%inibizione laterale
L = -1.2*ones(Nc,Nc)+diag(1.2*ones(Nc,1));   %tolto autoanello 1.2
%pesi da talamo a corteccia (sinapsi eccitatorie)
Wct = 4*diag(ones(Nc,1));  %diagonale

% %%%%%%%%%%%% addestrabili %%%%%%%%%%%%%%%%%%%%%
% 
% %inizializzazione
% punizione = 1.;
% caso = 0.;
% Nc = 4;
% Wgc = 0.24*diag(ones(Nc,1)+caso*randn(Nc,1))/punizione;
% Wgs = 0.9*diag(ones(Nc,1)+caso*randn(Nc,1))/punizione;
% Wnc = 1.08*diag(ones(Nc,1)+caso*randn(Nc,1))*punizione;
% Wns = 0.1*diag(ones(Nc,1)+caso*randn(Nc,1))*punizione;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pesi da No-Go a Gpe (sinapsi inibitorie)
Wen = -2*diag(ones(Nc,1));   %diagonale
Wen = Wen+10/100*Wen;   %diagonale     %%%%MODIFICA 22/12/2014

%pesi da Gpe a Gpi (sinapsi inibitorie)
Wie = -3*diag(ones(Nc,1));   %diagonale
%pesi da Go a Gpi (sinapsi inibitorie)
Wig = -36*diag(ones(Nc,1));   %diagonale  -12

%pesi da corteccia a talamo (sinapsi eccitatorie)
Wtc = 3*diag(ones(Nc,1));   %diagonale
%pesi da Gpi a talamo (sinapsi inibitorie)
Wti = -3*diag(ones(Nc,1));   %diagonale


%%%%%%%%%%%%%%%%%%% pesi da-a STN %%%%%%%%%%%%%%%%%%%
%peso da energia corteccia a STN (eccitazione)
%Ke = 7;    adesso � passato come parametro
%peso da Gpe a STN (inibizione)
Kgpe = -1;

%pesi da STN a Gpe (sinapsi eccitatorie)
Westn = 1;

%peso da STN a Gpi (sinapsi eccitatorie)
Wistn = 30;   %14;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% pesi da ChI a Go-NoGo %%%%%%%%%%%%%%%%%%%
%peso da interneurone colinergico a Go (sinapsi inibitorie)
wgchi = -1;

%pesi da interneurone colinergico a No-Go (sinapsi eccitatorie)
wnchi = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% guadagni %%%%%%%%%%%%%%%%%%%
%guadagno da DA a Go (eccitazione)
Ugo_trigger = 1.074;
% alpha = 0.75;  %(0.2*(Ugo_trigger-0.8)+0.5)/(0.7*(Ugo_trigger-0.8));
% 
% %guadagno da DA a No-Go (inibizione)
% beta = -1;
% 
% %guadagno da DA a interneurone colinergico (inibizione)
% gamma = -0.5;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calcolo attivit� dinamica delle strutture: dinamica (passabasso) + sigmoide

%parametri sigmoide
a = 4;
U0 = 1.0;

%attivit� tonica di Gpi > attivit� tonica di Gpe (ref.)
%attivit� tonica di Gpe
Igpe = 1.0;
%attivit� tonica di Gpi
Igpi = 3;

%attivit� tonica di ChI
Ichi = 1.00;

% assegno dopamina tonica+dopamina ext
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dop_tonic = 0.40;
DA = Dop_tonic;



%% condizioni iniziali


Ugpe(:,1) = Igpe;
Ugpi(:,1) = Igpi;
Uchi(1) = Ichi+gamma*DA;

% STN_ON = 1;
% T_ON = 1;

Ns = 4;
% 
% S = zeros(Ns,1);
% S(1) = 1;

t_alto = 0;
t_reset = 0;
vincitore = 0;
T_reset = 55;  %65;   %70
T_alto = 35;  %40;    %45
k_tap_vett = [];
caso_neuron = 0.0;
%%
S0(1) = 1;
S0(2:Nc) = 0;
noise1 =   caso_neuron*randn(Nc,1);
for k = 1:D
    
   
   
    if (t(k) > t_alto) && (vincitore == 1)   || ( (t(k) > t_reset + 500) && (vincitore == 0) ) %per 200 ms non ho risposto
    S0(1:2) = S0(1:2)*(-1) +1;   % cambio segno agli ingressi, ma solo dopo l'istante t_alto; S(1:2) = S(1:2)*(-1) +1
    S = S0';
    vincitore = 0;
    t_reset = t(k) + T_reset;
    end
    
    if t(k) < t_reset
        C(:,k) = 0;
        Uc(:,k) = 0;
%         Ugo(:,k) = 0;
%         Unogo(:,k) = 0;
%         Ut(:,k) = 0;
    else
        C(:,k) = 1./(1+exp(-a*(Uc(:,k)-U0)));
    end
    
    if (vincitore == 0) && ( max(C(:,k)) > 0.9)
        t_alto = t(k) + T_alto;
        k_tap_vett = [k_tap_vett t(k)];
        vincitore = 1;
    end
    
 
    
    %C(:,k) = 1./(1+exp(-a*(Uc(:,k)-U0)));
    
    Go(:,k) = 1./(1+exp(-a*(Ugo(:,k)-U0)));
    NoGo(:,k) = 1./(1+exp(-a*(Unogo(:,k)-U0)));
    Gpe(:,k) = 1./(1+exp(-a*(Ugpe(:,k)-U0)));
    Gpi(:,k) = 1./(1+exp(-a*(Ugpi(:,k)-U0)));
    ChI(k) = 1./(1+exp(-a*(Uchi(k)-U0)));
    
    if T_ON==1
        T(:,k) = 1./(1+exp(-a*(Ut(:,k)-U0)));
        %altrimenti resta a 0 come inizializzato
    elseif T_ON~=0
        disp('Wrong value for T_ON')
        return
    end
    
    if STN_ON==1
        STN(k) = 1./(1+exp(-a*(Ustn(k)-U0)));
        %altrimenti resta a 0 come inizializzato
    elseif STN_ON~=0
        disp('Wrong value for STN_ON')
        return
    end
    
    
   
    %energia in corteccia (per conflitto)
    for i = 1:Nc
        for j = i:Nc
            E(k) = E(k)+C(i,k)*C(j,k);
        end
    end
    E(k) = E(k)-(sum(C(:,k).^2));
    
    %%%% dinamica
    Ul(:,k+1) = Ul(:,k)+dt/tauL*(-Ul(:,k)+L*C(:,k));
    Uc(:,k+1) = Uc(:,k)+dt/tau*(-Uc(:,k)+Wcs*S+Ul(:,k)+Wct*T(:,k))+noise1;
    
    IGo_DA_Ach(:,k) = alpha*DA*(Go(:,k)-0.35)+wgchi*ChI(k); %+0.18; 
    INoGo_DA_Ach(:,k) = (beta*DA+wnchi*ChI(k))*ones(Nc,1);
    
    Ugo(:,k+1) = Ugo(:,k)+dt/tau*(-Ugo(:,k)+Wgs*S+Wgc*C(:,k)+IGo_DA_Ach(:,k));     %versione originaria
    Unogo(:,k+1) = Unogo(:,k)+dt/tau*(-Unogo(:,k)+Wns*S+Wnc*C(:,k)+INoGo_DA_Ach(:,k));
    Ugpe(:,k+1) = Ugpe(:,k)+dt/tau*(-Ugpe(:,k)+Wen*NoGo(:,k)+Westn*STN(k)+Igpe);
    Ugpi(:,k+1) = Ugpi(:,k)+dt/tau*(-Ugpi(:,k)+Wig*Go(:,k)+Wie*Gpe(:,k)+Wistn*STN(k)+Igpi);
    Uchi(k+1) = Uchi(k)+dt/tau*(-Uchi(k)+Ichi+gamma*DA);   %+gammaDop_tonic
    
    
    if T_ON==1
        Ut(:,k+1) = Ut(:,k)+dt/tau*(-Ut(:,k)+Wti*Gpi(:,k)+Wtc*C(:,k));
        %altrimenti resta a 0 come inizializzato
    elseif T_ON~=0
        disp('Wrong value for STN_ON')
        return
    end
    
    if STN_ON==1
        Ustn(k+1) = Ustn(k)+dt/tau*(-Ustn(k)+Ke*E(k)+Kgpe*sum(Gpe(:,k)));
        %altrimenti resta a 0 come inizializzato
    elseif STN_ON~=0
        disp('Wrong value for STN_ON')
        return
    end
    %%%%
    
end


%% calcolo frequenza tapping

if length(k_tap_vett)< 3
    ft = 0;
else
Tt = k_tap_vett(end) - k_tap_vett(end-2);
ft = 1/Tt*1000;
end
frequenza = ft*60;





