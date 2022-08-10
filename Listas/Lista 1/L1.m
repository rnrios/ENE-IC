clear 
close all
clc

%%% Questão 1 %%%

% Mi = 227;
% Sigma = 3;
% x = Sigma.*randn(100,1) + Mi;
% stats = [mean(x) std(x)]
% nbins = 10;
% h = histogram(x,nbins)
% title("Figura 1")


%%% Questão 2 %%%

% x = unifrnd(0,10,[1000 1]);
% h1 = histogram(x,12);
% title("Figura 2 b)")
% x_t = reshape(x,[50 20])';
% x_a = zeros(20,1);
% for i=1:20
%     x_a(i) = mean(x_t(i,:));
% end
% figure
% min(x_a)
% max(x_a)
% h2 = histogram(x_a,6);
% title("Figura 2 c)")


%%% Questão 5 %%%

Kb = 3.1e6; Bb = -2254.473; dK = 1.2e5; dB = 11.94;
Omax = 1e4; Omin = 1e2; Imax = 393; Imin = 218;

% Item a
k = (Omax-Omin)/(Imax-Imin);
a = Omin - k*Imin;
% Item b
syms T
Rideal = k*T+a;
R = Kb*exp(Bb/T);
N = vpa(R - Rideal,4);
dN = vpa(diff(N,T),4);
d2N = vpa(diff(dN,T),4);
I = Imin:0.1:Imax;
d2N = subs(d2N,T,I);
% plot(I,d2N)
% title("Segunda Derivada de N(T)")
%Item c
T_c = solve(dN == 0,T,'Real',true);
T_c = T_c(1);
N_c = subs(N,T,T_c); 
N_c/(Omax-Omin)*100;
% Item d
vpa(R,4);
dR = vpa(diff(R,T),4);
subs(dR,T,T_c);
%Item e
syms K B
Rs = double(subs(R,T,I));
R = K*exp(B/T);
sk = 1.2e5; sb = 11.94;
dRB = vpa(subs(diff(R,B),[B K T],[Bb Kb 373]),4);
dRK = vpa(subs(diff(R,K),[B K T],[Bb Kb 373]),4);
dR = sqrt(dRB^2*sb^2+dRK^2*sk^2);
% %Item f
Rideals = subs(Rideal,T,I);
Ns = double(subs(N,T,I));
dNs = diff(Ns)./diff(I);
k = round((T_c-Imin)/0.1+1);
tang=(I-I(k))*dNs(k)+Ns(k);
plot(I,Rs,'LineWidth',2)
hold on
plot(I,Rideals,'LineWidth',2)
hold on
plot(I,Ns,'LineWidth',2,'Color',[0.4660 0.6740 0.1880])
hold on
plot(I,tang,'LineWidth',2)
legend("Termistor","Reta Ideal","Não-Linearidade","Tangente")
grid
title("Figura 5 f")


%%% Questão 6 %%%
% R0b = 100; ab= 3.91e-3; sR0 = 0.29; sa = 2.94e-5;
% syms Rt R0 a T
% %Item a
% Tb = 37;
% Rt = R0*(1+a*T);
% Rtb = vpa(subs(Rt,[R0 a T],[R0b ab Tb]),4);
% dR0 = vpa(subs(diff(Rt,R0),[R0 a T],[R0b ab Tb]),4);
% dRa = vpa(subs(diff(Rt,a),[R0 a T],[R0b ab Tb]),4);
% dRt = vpa(sqrt(dR0^2*sR0^2+dRa^2*sa^2),4);
% 
% syms Vs r Rt R1
% Vsb=10; R1b=100; rb=0.01; sVs=0; sR1 = 0.89; sr = 1.94e-4;
% V = Vs*r*(Rt/R1-1);
% Vb = Vsb*rb*(Rtb/R1b-1);
% dR1 = vpa(subs(diff(V,R1),[Vs R1 r Rt],[Vsb R1b rb Rtb]),4)
% dr = vpa(subs(diff(V,r),[Vs R1 r Rt],[Vsb R1b rb Rtb]),4)
% delRt = vpa(subs(diff(V,Rt),[Vs R1 r Rt],[Vsb R1b rb Rtb]),4)
% dV = sqrt(delRt^2*dRt^2+dR1^2*sR1^2+dr^2*sr^2)
% 
% syms n K1 V b1 
% K1b = 6522; b1b = 0; sk1 = 0; sb1 = 0.5/sqrt(3);
% sK1 = 0;
% n = K1*V + b1;
% nb = round(K1b*Vb + b1b);
% db1 = vpa(subs(diff(n,b1),[K1 V b1],[K1b Vb b1b]),4);
% delV = vpa(subs(diff(n,V),[K1 V b1],[K1b Vb b1b]),4);
% dn = sqrt(db1^2*sb1^2+dV^2*delV^2);
% 
% K2 = 0.391; b2 = 0;
% Tm = vpa(K2*nb,4);
% eTm = sqrt(K2^2*dn^2);


%%% Questão 7 %%%
% x = L1_import("Lista1_Q7.csv",1,1001);
% t = table2array(x(:,1));
% y = table2array(x(:,2));
% 
% Reg = stepinfo(y,t);
% ts = Reg.SettlingTime;
% tp = Reg.PeakTime
% Ov = Reg.Overshoot
% csi = -log(Ov/100)/sqrt(pi^2+(log(Ov/100))^2)
% wd = pi/tp
% wn = wd/sqrt(1-csi^2)
% G = tf(wn^2,[1 2*csi*wn wn^2])
% [y1,t1] = step(G,10);
% plot(t,y,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% plot(t1,squeeze(y1),'LineWidth',2)
% legend("Resposta Original","Resposta Item a")
% grid
% title("Figura 7b")


%%% Questão 8 %%%

% syms Ka
% T = 10; Tf = 2; K = 99;
% Ka = double(solve(Tf/(Ka+1)==0.2));
% G = tf(K,[T 1]);
% H1 = tf(1,[Tf 1]);
% H2 = tf(Ka,[Tf Ka+1]);
% Gmf1 = feedback(G,H1);
% Gmf2 = feedback(G,H2);
% stepinfo(Gmf1)
% stepinfo(Gmf2)


%%% Questão 9 %%%

% Rn = 1e5; Rc = 500; Rr = 250; Pin = 5e3;
% R = Rn/(Rn+Rc+Rr); Omax = 20e-3; Omin = 4e-3; Imax = 1e4; Imin = 0;
% k = (Omax-Omin)/(Imax-Imin)
% a = Omin - k*Imin
% i = k*Pin + a
% Vv = Rr*i
% Vin = R*Vv
% V = 1e4/4;
% Pb = -V;
% P = V*Vin + Pb
% E = P - Pin

%%% Questão 10 %%%
% Rin = 1e9; Rl = 1e4; Ro = 100;
% syms Rth
% Rth = solve(224*Rin/(Rth+Rin)*Rl/(Ro+Rl) == 5);
% i = 224/Rth;
% R1 = 5/i^2;
% syms R2
% R2 = solve((2*R1*R2)/(2*R1+R2) == Rth);
% vpa(i,4);
% vpa(R1,4);
% vpa(R2,4);
% vpa(R2*i^2,4);