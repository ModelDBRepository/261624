function costo = stima_levodopa(p)
global V1 V2 V3 t_dop y i ke3 k31 alpha beta gamma
global S Wgc Wgs Wnc Wns Ke STN_ON T_ON 
global c1 c3_delay tfreq Tapping_f Dop_tonic


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
%c1(1) = D/V1;


for j = 1: L-1,
    dc1 = -k21/V1*c1(j)+k12/V1*c2(j)-ke1/V1*c1(j)-k31/V1*c1(j)+i(j)/V1;
    dc2= k21/V2*c1(j)-k12/V2*c2(j);
    c1(j+1) = c1(j) +dt*dc1;
    c2(j+1) = c2(j) + dc2*dt;
end

    
    
% plot(t,y,'*-.r',t1,c2,'-b',t1,c3,'g')
% pause(0.1)
index = [1 151 301 451 601 751 901 1201 1501 1801];  
kk = length(y);
index = index(1:kk);  % inserito perché in qualche caso y è meno lungo
costo = sum((c1(index) - y).^2,'omitnan');
if min(p) < 0.1 costo = 10000000; end
if max(p) > 10 costo = 10000000; end






