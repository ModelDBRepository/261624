
function [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI] = BG_model_function_Ach(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Dop_tonic)
% BG_model_function_Ach                         returns dynamical behaviour of Basal Ganglia structures and cortex
% S                                             stimulus. Column vector of 4 elements Ei, 0<=Ei<=1 for i=1,2,3,4
% Wgc,Wgs,Wnc,Wns                               sets corresponding synaptic weights
% C_CORRECT_WINNER                              sets the desired correct response, if possible, depending on the stimulus
% Uc,Ugo,Unogo,Ugpe,Ugpi,Ut,Ustn                returns the input to the sigmoidal function within time of the corresponding brain structures
% C,Go,NoGo,Gpe,Gpi,T,STN,E                     returns activity within time of the corresponding brain structures and energy in the cortex
% IGo_DA_Ach,INoGo_DA_Ach                       returns the input due to Dopa and Ach to Go and NoGo units
% t                                             returns time
% Wgc_post,Wgs_post,Wnc_post,Wns_post           returns corresponding synaptic weights after Hebbian learning
% r                                             returns +1 for reward, -1 for punishment, NaN for no feedback
% r_story                                       retyrn the historical record of rewards and punishments
% k_reward                                      returns position of feedback, NaN for no feedback
% ChI                                           returns activity within time of the cholinegic interneuron

%% tempi

tau = 10;%10;   %costante di tempo di base
tauL = 5*tau;   %costante di tempo dell'inibizione laterale
%tauS = tau/5;

dt = 0.1;   %passo
t = (0:dt:800)';   %tempi [ms]
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

%pesi da stimolo a corteccia
Wcs = 1*ones(4,4);   %0.2 extradiag, 1.1 su diag
%inibizione laterale
L = -1.2*ones(Nc,Nc)+diag(1.2*ones(Nc,1));   %tolto autoanello
%pesi da talamo a corteccia (sinapsi eccitatorie)
Wct = 4*diag(ones(Nc,1));  %diagonale

%%%%%%%%%%%% addestrati %%%%%%%%%%%%%%%%%%%%%

% %pesi da corteccia a Go (sinapsi eccitatorie)
% %Wgc = 0.4*diag(ones(Nc,1));   %diagonale
% Wgc = 0.6*diag(ones(Nc,1));
% %pesi da stimolo a Go
% %Wgs = 0.4*diag(ones(Nc,1));
% Wgs = 0.6*diag(ones(Nc,1));
%
% %pesi da corteccia a No-Go (sinapsi eccitatorie)
% Wnc = 1*diag(ones(Nc,1));
% %Wnc = 1.2*diag(ones(Nc,1));
% %pesi da stimolo a No-Go
% Wns = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pesi da No-Go a Gpe (sinapsi inibitorie)
Wen = -2.2*diag(ones(Nc,1));   %diagonale

%pesi da Gpe a Gpi (sinapsi inibitorie)
Wie = -3*diag(ones(Nc,1));   %diagonale
%pesi da Go a Gpi (sinapsi inibitorie)
Wig = -36*diag(ones(Nc,1));   %diagonale

%pesi da corteccia a talamo (sinapsi eccitatorie)
Wtc = 3*diag(ones(Nc,1));   %diagonale
%pesi da Gpi a talamo (sinapsi inibitorie)
Wti = -3*diag(ones(Nc,1));   %diagonale

%%%%%%%%%%%%%%%%%%% pesi da-a STN %%%%%%%%%%%%%%%%%%%
%peso da energia corteccia a STN (eccitazione)
Ke = 7;
%peso da Gpe a STN (inibizione)
Kgpe = -1;

%pesi da STN a Gpe (sinapsi eccitatorie)
Westn = 1;

%peso da STN a Gpi (sinapsi eccitatorie)
Wistn = 30;   %14
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
alpha = 0.75;  %11; %(0.2*(Ugo_trigger-0.8)+0.5)/(0.7*(Ugo_trigger-0.8));

%guadagno da DA a No-Go (inibizione)
beta = -1;

%guadagno da DA a interneurone colinergico (inibizione)
gamma = -0.5;   %-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% calcolo attività dinamica delle strutture: dinamica (passabasso) + sigmoide

%parametri sigmoide
a = 4;
U0 = 1.0;

%attività tonica di Gpi > attività tonica di Gpe (ref.)
%attività tonica di Gpe
Igpe = 1.0;
%attività tonica di Gpi
Igpi = 3;

%attività tonica di ChI
Ichi = 1.00;   % 1.25

%dopamina tonica in sano
% Dop_tonic_CTR = 0.55;    %0.55
% Dop_tonic = Dop_tonic_CTR;

%reward o punishment
%r = NaN   no reward o inizializzazione
%r = 1   reward
%r = -1   punishment
r = NaN;
%istante di reward o punishment
%k_reward = NaN   no reward o inizializzazione
k_reward = NaN;


%tempistiche dell'erogazione dopamina fasica
latency = 100;  %ms, valore fisiologico 50-110 ms (Schultz 1998)
klatency = ceil(latency/dt);
duration = 50;   %ms, valore fisiologico <200 ms (Schultz 1998)
kduration = ceil(duration/dt);


%% condizioni iniziali


Ugpe(:,1) = Igpe;
Ugpi(:,1) = Igpi;
Uchi(1) = Ichi+gamma*Dop_tonic;

Ns = 4;

r_story = [];
noise1 =   0.15*randn(Nc,1);

for k = 1:D
    
    %EROGAZIONE DOPAMINA FASICA
    
    if isnan(r);   %nulla
        dop_phasic = 0;
    else   %k_reward è stato inizializzato
                if k>=k_reward+klatency && k<=k_reward+klatency+kduration   %sono in arco di tempo di variazione dopamina
                    if r==1   %reward
                        delta_Dop = 1.2*Dop_tonic;
                    elseif r==-1   %punishment
                        delta_Dop = -Dop_tonic;
                    end
                    dop_phasic = delta_Dop;
                else   %sono fuori arco di tempo di variazione dopamina
                    dop_phasic = 0;
                end
    end

    DA = Dop_tonic+dop_phasic;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% sigmoide
    C(:,k) = 1./(1+exp(-a*(Uc(:,k)-U0)));
    Go(:,k) = 1./(1+exp(-a*(Ugo(:,k)-U0)));
    NoGo(:,k) = 1./(1+exp(-a*(Unogo(:,k)-U0)));
    Gpe(:,k) = 1./(1+exp(-a*(Ugpe(:,k)-U0)));
    Gpi(:,k) = 1./(1+exp(-a*(Ugpi(:,k)-U0)));
    ChI(k) = 1./(1+exp(-a*(Uchi(k)-U0)));
    T(:,k) = 1./(1+exp(-a*(Ut(:,k)-U0)));
    STN(k) = 1./(1+exp(-a*(Ustn(k)-U0)));
    %%%%
    
        %REWARD, PUNISHMENT O NIENTE
        %vedo se c'è 1 vincitore, re c'è devo erogare reward o punishment
        %lo faccio solo una volta, basta che raggiunga vincita, per evitare
        %oscillazioni poi dovute all'erogazione del reward o punishment
        if isnan(Correct_winner)==0;   %ho modo di dare feedback, so chi deve vincere  %MODIFICA VERSIONE 9!!!!
            if isnan(r);   %non ho ancora un esito della prova, quindi procedo a verificare
                C_winner_list = find(C(:,k)>=0.9);
                lost = isempty(C_winner_list);
                if lost==0   %c'è almeno un vincitore
                    n_C_winner = length(C_winner_list);
                    if n_C_winner==1   %c'è un solo vincitore, procedo al reward o punishment
                        k_reward = k;   %prendo istante in cui si verifica la vittoria
                        if C_winner_list==Correct_winner
                            r = 1;   %reward
                        else
                            r = -1;   %punishment
                        end
                        r_story = [r_story; r];
                    end
                end
            end
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
    Uc(:,k+1) = Uc(:,k)+dt/tau*(-Uc(:,k)+Wcs*S+Ul(:,k)+Wct*T(:,k)+noise1);
    
    IGo_DA_Ach(:,k) = alpha*DA*(Go(:,k)-0.35)+wgchi*ChI(k);
    INoGo_DA_Ach(:,k) = (beta*DA+wnchi*ChI(k))*ones(Nc,1);
    
    Ugo(:,k+1) = Ugo(:,k)+dt/tau*(-Ugo(:,k)+Wgs*S+Wgc*C(:,k)+IGo_DA_Ach(:,k));     %versione originaria
    
    %Ugo(:,k+1) = (Wgs*S+Wgc*C(:,k)-0.7*alpha*0.8+wgchi*ChI(k))/(1-alpha*0.7);     %riscrtitta il 09/12/2014
    
    Unogo(:,k+1) = Unogo(:,k)+dt/tau*(-Unogo(:,k)+Wns*S+Wnc*C(:,k)+INoGo_DA_Ach(:,k));
%     Ugo(:,k+1) = Ugo(:,k)+dt/tau*(-Ugo(:,k)+Wgs*S+Wgc*C(:,k)+alpha*DA*(Ugo(:,k)-0.8)+wgchi*ChI(k));     
%     Unogo(:,k+1) = Unogo(:,k)+dt/tau*(-Unogo(:,k)+Wns*S+Wnc*C(:,k)+beta*DA+wnchi*ChI(k));
    Ugpe(:,k+1) = Ugpe(:,k)+dt/tau*(-Ugpe(:,k)+Wen*NoGo(:,k)+Westn*STN(k)+Igpe);
    Ugpi(:,k+1) = Ugpi(:,k)+dt/tau*(-Ugpi(:,k)+Wig*Go(:,k)+Wie*Gpe(:,k)+Wistn*STN(k)+Igpi);
    Uchi(k+1) = Uchi(k)+dt/tau*(-Uchi(k)+Ichi+gamma*DA);
    Ut(:,k+1) = Ut(:,k)+dt/tau*(-Ut(:,k)+Wti*Gpi(:,k)+Wtc*C(:,k));
    Ustn(k+1) = Ustn(k)+dt/tau*(-Ustn(k)+Ke*E(k)+Kgpe*sum(Gpe(:,k)));
    %%%%
    
end

%% REGOLA DI HEBB
%inizializzazione
S_th = 0.5*ones(length(S),1);
C_th = 0.5*ones(Nc,1);
Go_th = 0.5*ones(Nc,1);
NoGo_th = 0.5*ones(Nc,1);

%gamma = 0.05;
%gain = 0.1;  %OK X RISULTATI IN MEMORIA
%gamma = 0.001;
gain = 0.01;


delta_Wgc = zeros(Nc,Nc);
delta_Wgs = zeros(Nc,Nc);

delta_Wnc = zeros(Nc,Nc);
delta_Wns = zeros(Nc,Nc);


if isnan(r)==0   %c'è stato reward o punishment
    
    if ((k_reward+klatency)>D)==0 &&((k_reward+klatency+kduration)>D)==0   %non addestro se andamento sta variando o sta finendo di variare
%         if r==1   %reward
%             %massimo del Go che ha vinto e minimo dei Go che hanno perso
%             %minimo di tutti i NoGo
%             val_Go = min(Go(:,(k_reward+klatency):(k_reward+klatency+kduration)),[],2);
%             val_Go(C_winner_list) = max(Go(C_winner_list,(k_reward+klatency):(k_reward+klatency+kduration)),[],2);
%             val_NoGo = min(NoGo(:,(k_reward+klatency):(k_reward+klatency+kduration)),[],2);
%         else   %punishment
%             %minimo di tutti i Go
%             %massimo di tutti i NoGo
%             val_Go = min(Go(:,(k_reward+klatency):(k_reward+klatency+kduration)),[],2);
%             val_NoGo = max(NoGo(:,(k_reward+klatency):(k_reward+klatency+kduration)),[],2);
%         end

%% cambio
if r==1   %reward
    %massimo del Go che ha vinto e minimo dei Go che hanno perso
    %minimo di tutti i NoGo
    val_Go(1) = min(Go(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(2) = min(Go(2,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(3) = min(Go(3,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(4) = min(Go(4,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(C_winner_list) = max(Go(C_winner_list,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(1) = min(NoGo(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(2) = min(NoGo(2,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(3) = min(NoGo(3,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(4) = min(NoGo(4,(k_reward+klatency):(k_reward+klatency+kduration)));
else   %punishment
    %minimo di tutti i Go
    %massimo di tutti i NoGo
    val_Go(1) = min(Go(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(2) = min(Go(2,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(3) = min(Go(3,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_Go(4) = min(Go(4,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(1) = max(NoGo(1,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(2) = max(NoGo(2,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(3) = max(NoGo(3,(k_reward+klatency):(k_reward+klatency+kduration)));
    val_NoGo(4) = max(NoGo(4,(k_reward+klatency):(k_reward+klatency+kduration)));
end


        
%%
        for l=1:Nc
            delta_Wgc(l,l) = gain*max(0,(C(l,k_reward)-C_th(l))).*(val_Go(l)-Go_th(l));
            delta_Wnc(l,l) = gain*max(0,(C(l,k_reward)-C_th(l))).*(val_NoGo(l)-NoGo_th(l));
        end
        
        for row=1:Nc
            for col=1:Nc
                delta_Wgs(row,col) = gain*(val_Go(row)-Go_th(row))*max(0,(S(col)-S_th(col)));
                delta_Wns(row,col) = gain*(val_NoGo(row)-NoGo_th(row))*max(0,(S(col)-S_th(col)));
            end
        end
    end
    
end

Wgc_post = Wgc+delta_Wgc;
Wgs_post = Wgs+delta_Wgs;
Wnc_post = Wnc+delta_Wnc;
Wns_post = Wns+delta_Wns;


%% CONTROLLI SULLE SINAPSI

%saturazione sinapsi

%sinapsi eccitatorie

%delta_opt = 0.05;
%delta_opt = 0.15;
delta_opt = 0.2;

sinapsi_ecc_max = 1.2;
max_ecc = sinapsi_ecc_max+delta_opt;

Wgc_opt_max = max_ecc*diag(ones(Nc,1));   %diagonale
Wgc_opt_min = zeros(Nc,Nc);
Wnc_opt_max = max_ecc*diag(ones(Nc,1));   %diagonale
Wnc_opt_min = zeros(Nc,Nc);
Wgs_opt_max = max_ecc*ones(Nc,Nc);   %piena
Wgs_opt_min = zeros(Nc,Nc);
Wns_opt_max = max_ecc*ones(Nc,Nc);   %piena
Wns_opt_min = zeros(Nc,Nc);

[rowgc,colgc,valgc] = find(Wgc_post>Wgc_opt_max);
if isempty(rowgc)==0
    for i=1:length(rowgc)
        Wgc_post(rowgc(i),colgc(i)) = Wgc_opt_max(rowgc(i),colgc(i));
    end
end
[rowgc,colgc,valgc] = find(Wgc_post<Wgc_opt_min);
if isempty(rowgc)==0
    for i=1:length(rowgc)
        Wgc_post(rowgc(i),colgc(i)) = Wgc_opt_min(rowgc(i),colgc(i));
    end
end

[rownc,colnc,valnc] = find(Wnc_post>Wnc_opt_max);
if isempty(rownc)==0
    for i=1:length(rownc)
        Wnc_post(rownc(i),colnc(i)) = Wnc_opt_max(rownc(i),colnc(i));
    end
end
[rownc,colnc,valgnc] = find(Wnc_post<Wnc_opt_min);
if isempty(rownc)==0
    for i=1:length(rownc)
        Wnc_post(rownc(i),colnc(i)) = Wnc_opt_min(rownc(i),colnc(i));
    end
end

[rowgs,colgs,valgs] = find(Wgs_post>Wgs_opt_max);
if isempty(rowgs)==0
    for i=1:length(rowgs)
        Wgs_post(rowgs(i),colgs(i)) = Wgs_opt_max(rowgs(i),colgs(i));
    end
end
[rowgs,colgs,valgs] = find(Wgs_post<Wgs_opt_min);
if isempty(rowgs)==0
    for i=1:length(rowgs)
        Wgs_post(rowgs(i),colgs(i)) = Wgs_opt_min(rowgs(i),colgs(i));
    end
end

[rowns,colns,valns] = find(Wns_post>Wns_opt_max);
if isempty(rowns)==0
    for i=1:length(rowns)
        Wns_post(rowns(i),colns(i)) = Wns_opt_max(rowns(i),colns(i));
    end
end
[rowns,colns,valns] = find(Wns_post<Wns_opt_min);
if isempty(rowns)==0
    for i=1:length(rowns)
        Wns_post(rowns(i),colns(i)) = Wns_opt_min(rowns(i),colns(i));
    end
end
