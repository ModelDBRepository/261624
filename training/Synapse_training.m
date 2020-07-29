

%% ADDESTRAMENTO

%% stimoli base

%Sn: stimolo, ipotesi di diversi stimoli applicati alla rete
Ns = 4;

%S1: stimolo 1
S1 = zeros(Ns,1);
S1(1) = 1.0;   %stimolo massimo in posizione 1
S1(2) = 0.;   %stimolo alto in posizione 2
S1(3) = 0;
S1(4) = 0;

Correct_winner_1 = 1;

%S2: stimolo 2
S2 = zeros(Ns,1);
S2(1) =  0.;   %stimolo alto in posizione 1
S2(2) = 1.0;   %stimolo massimo in posizione 2
S2(3) = 0;
S2(4) = 0;

Correct_winner_2 = 2;


%% %%%%%%%%%%%%%%%%%% SOGGETTO NORMALE %%%%%%%%%%%%%%%%%%%%%%%%%%

%% inizializzazione sinapsi

Nc = 4;
%     
%     
%     par = 1/10;
%     %pesi da corteccia a Go (sinapsi eccitatorie)
%     Wgc = 0.48*diag(ones(Nc,1))-par*0.48*diag(ones(Nc,1));
%     Wgc(3,3) = 0;
%     Wgc(4,4) = 0;
%     
%     
%     %pesi da stimolo a Go (sinapsi eccitatorie)
%     Wgs = 0.6*diag(ones(Nc,1))-par*0.6*diag(ones(Nc,1));
%     Wgs(1,2) = 0.2;
%     Wgs(2,1) = 0.2;
%     
%     %pesi da corteccia a NoGo (sinapsi eccitatorie)
%     Wnc = 1.08*diag(ones(Nc,1))-par*1.08*diag(ones(Nc,1));
%     
%     %pesi da stimolo a NoGo (sinapsi eccitatorie)
%     Wns = 0.4*diag(ones(Nc,1))-par*0.4*diag(ones(Nc,1));
%     Wns(1,2) = 0.2;
%     Wns(2,1) = 0.2;

%%

N_epochs = 50;  %PER NaN

j1_reward = 3;
j1_punishment = 5;
% S1_reward = 3:4:(3+4*100);
% S1_punishment = 5:4:(5+4*100);

j2_reward = 2;
j2_punishment = 4;
% S2_reward = 2:4:(2+4*100);
% S2_punishment = 4:4:(2+4*100);

Wgc_epocs = zeros(Nc,Nc,N_epochs);
Wgs_epocs= zeros(Nc,Nc,N_epochs);
Wnc_epocs = zeros(Nc,Nc,N_epochs);
Wns_epocs = zeros(Nc,Nc,N_epochs);

Wgc_epocs(:,:,1) = Wgc;
Wgs_epocs(:,:,1) = Wgs;
Wnc_epocs(:,:,1) = Wnc;
Wns_epocs(:,:,1) = Wns;

vett_reward = zeros(N_epochs,1);
vett_punishment = zeros(N_epochs,1);
vett_no_risposta = zeros(N_epochs,1);

S_vett = zeros(2,N_epochs);
S1_vett = zeros(2,N_epochs);
S2_vett = zeros(2,N_epochs);



%%
for i = 1:N_epochs
    
    i
    Wgc = squeeze(Wgc_epocs(:,:,i));
    Wgs = squeeze(Wgs_epocs(:,:,i));
    Wnc = squeeze(Wnc_epocs(:,:,i));
    Wns = squeeze(Wns_epocs(:,:,i));


    resto = rem(i,2);
    noise = 0.*randn(2,1);

    if resto == 1   %sono nel caso dispari
        S = S1;  
        S(1) = S(1)+noise(1);
        S(2) = S(2)+noise(2);
        Correct_winner = Correct_winner_1;
    elseif resto == 0   %sono nel caso pari
        S = S2;
        S(1) = S(1)+noise(1);
        S(2) = S(2)+noise(2);
        Correct_winner = Correct_winner_2;
    end
    
    S(find(S>1)) = 1;
    S(find(S<0)) = 0;
    S(3:4) = 0;
    
    
    
    [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,tt,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI] = BG_model_function_Ach(S,Wgc,Wgs,Wnc,Wns,Correct_winner,DA);


%     i 
%     r
    
    if r==1
        vett_reward(i) = 1;
    elseif r==-1
        vett_punishment(i) = 1;
    else
        vett_no_risposta(i) = 1;
    end
    
    S_vett(1,i) = S(1);
    S_vett(2,i) = S(2);
    
    if resto == 1
        S1_vett(1,i) = S(1);
        S1_vett(2,i) = S(2);

    elseif resto == 2
        S2_vett(1,i) = S(1);
        S2_vett(2,i) = S(2);
    end

    
    Wgc_epocs(:,:,i+1) = Wgc_post;
    Wgs_epocs(:,:,i+1) = Wgs_post;
    Wnc_epocs(:,:,i+1) = Wnc_post;
    Wns_epocs(:,:,i+1) = Wns_post;
    
    

    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    clear Wgc Wgs Wnc Wns
    
end 

%premi
reward_tot = sum(vett_reward)
%punizioni
punishment_tot = sum(vett_punishment)
%no risposte
no_answer_tot = sum(vett_no_risposta)




