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
    dc3= k31/V3*c1(j) -ke3/V3*c3(j);
    c3(j+1) = c3(j) + dc3*dt;
end






