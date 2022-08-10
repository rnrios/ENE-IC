clear
close all
clc

%%%%%%%% Questão 2.1 %%%%%%%%
% syms a1 a2 a3
% T = [100 419.527 961.78];
% E = [645; 3375; 9149];
% 
% A = [T(1) T(1)^2 T(1)^3;T(2) T(2)^2 T(2)^3;T(3) T(3)^2 T(3)^3];
% cts = A\E;


%%%%%%%% Questão 2.2 %%%%%%%%
% syms a b
% 
% R = [9 .5];
% Theta = [273.15 373.15];
% 
% eqn = [R(1) == a*exp(b/Theta(1)),R(2) == a*exp(b/Theta(2))];
% vars = [a,b];
% [a,b] = solve(eqn,vars);
% vpa(a,4)
% vpa(b,4)
% vpa(a*exp(b/(298.15)),4)

%%%%%%%% Questão 2.6 %%%%%%%%

% x = [0 100 419.6];
% y = [100; 138.5; 253.7];
% A = [1 x(1) x(1)^2;1 x(2) x(2)^2;1 x(3) x(3)^2];
% p = vpa(A\y,4)


%%%%%%%% Questão 3.1 %%%%%%%%
% E = [0 4.017e-2 4.66e-6]*[1 117 117^2]';
% 
% K1 = 3.893;
% DTa = -10;
% a1 = -3.864;
% Km = 1.95e-4;
% Ki = 2e-3;
% i = K1*E+Km*E*DTa+Ki*DTa+a1;
% 
% K2 = 6.25;
% a2 = 25;
% Tb = K2*i+a2;
% E = Tb - 117;
% 
% dE = 6.93e-2;
% di = sqrt((K1+Km*DTa)^2*dE^2+0.14^2+(Km*E+Ki)^2*100);
% sqrt(K2^2*di^2+0.3^2)

%%%%%%%% Questão 3.2 %%%%%%%%

% dI2 = sqrt(0.005^2/3);
% dI3 = sqrt(4e-2^2*dI2^2+ 0.0005^2/3);
% dI4 = sqrt(1e3^2*dI3^2+ 0.5^2/3);
% dI5 = sqrt(225^2*dI4^2+ 100^2/3)

%%%%%%%% Questão 3.3 %%%%%%%%
% Vs = 1.5;
% F = 50;
% 
% V0 = 10*Vs/(1+100*Vs)*F;

%%%%%%%% Questão 3.4 %%%%%%%%

% syms Ks
% 
% Theta = 100*Ks/(1+10000*Ks);
% x = vpa(subs(Theta,5e-2),4);
% y = vpa(subs(Theta,1.1*5e-2),4);
% y-x

%%%%%%%% Questão 3.5 %%%%%%%%
% a1 = 4.3796e-2;
% a2 = -1.7963e-5;
% T1 = 100;
% T2 = 20;
% E = a1*(T1-T2)+a2*(T1^2-T2^2);
% 
% K1 = 255;
% V = K1*E;
% 
% K2 = 0.1;
% n = K2*V;
% 
% K3 = 1;
% b3 = 20;
% Tm = K3*n + b3;
% Erro = Tm -100
% 
% dE = sqrt(0.05^2 +(-a1-2*a2*T2)^2*2^2);
% dV = sqrt(5^2 + K1^2*dE^2);
% dn = sqrt(K2^2*dV^2 + (0.5/sqrt(3))^2);
% dTm = dn

%%%%%%%% Questão 3.6 %%%%%%%%

% r = 1.2;
% vt = 14;
% P = 0.5*r*vt^2;
% 
% K1 = 0.064;
% a1 = 4;
% i = K1*P + a1;
% 
% K2 = 12.8;
% n = round(K2*i);
% 
% K3 = 1.43;
% Vm = K3*sqrt(n-51);
% E = Vm - 14
% 
% dP = sqrt((0.5*vt^2)^2*(0.1/sqrt(3))^2);
% di = sqrt(K1^2*dP^2 + (0.04/sqrt(3))^2);
% dn = sqrt(K2^2*di^2 + (0.5/sqrt(3))^2);
% dVm = sqrt((K3/(2*sqrt(n-51)))^2*dn^2)


%%%%%%%% Questão 3.7 %%%%%%%%

% % Item a)
% Te = 320;
% K1m = 5e-4;
% Bm = 3e3;
% Vsm = -3;
% a1m = 0.77;
% K2m = 50;
% a2m = 300;
% 
% Rm = K1m*exp(Bm/Te);
% V0m = Vsm*(1/(1+3.3/Rm)-a1m);
% Tmm = K2m*V0m + a2m;
% E = Tmm - 320;
% 
% % Item b)
% syms T K1 B Vs a1 K2 a2 R V0 Tm
% 
% R = K1*exp(B/T);
% delR = diff(R,K1);
% delR = subs(delR,{B,T},{Bm,Te});
% dR = sqrt(delR^2*0.5e-4^2);
% 
% syms R
% V0 = Vs*(1/(1+3.3/R)-a1);
% delV0_1 = diff(V0,Vs);
% delV0_1 = subs(delV0_1,{Vs,R,a1},{Vsm,Rm,a1m});
% delV0_2 = diff(V0,R);
% delV0_2 = subs(delV0_2,{Vs,R,a1},{Vsm,Rm,a1m});
% delV0_3 = diff(V0,a1);
% delV0_3 = subs(delV0_3,{Vs,R,a1},{Vsm,Rm,a1m});
% dV0 = sqrt(delV0_1^2*0.03^2 + delV0_2^2*dR^2 + delV0_3^3*0.01^2);
% 
% Tm = K2*V0 + a2;
% delTm = sqrt(K2m^2*dV0^2 + 9);
% vpa(delTm,4)


%%%%%%%% Questão 3.8 %%%%%%%%

% -(1-1e-2*5e-2*1e3*1.9)*10;

%%%%%%%% Questão 3.9 %%%%%%%%

% (0.05*21.5*0.99 - 1)*5