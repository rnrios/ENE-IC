clear;
close all;
clc;

%%% Q1 %%%
% syms a1 a2 a3 T T2
% E = a1*T+a2*T^2+a3*T^3-a1*T2;
% Emin = vpa(subs(E,[a1 a2 a3 T2 T],[50.37 3.043e-2 -8.567e-5 0 -20]),4);
% Emax = vpa(subs(E,[a1 a2 a3 T2 T],[50.37 3.043e-2 -8.567e-5 0 100]),4);
% k = (Emax-Emin)/120;
% a = Emax-(100)*k;
% Ei = k*T + a;
% E_T = vpa(subs(E,[a1 a2 a3 T2],[50.37 3.043e-2 -8.567e-5 0]),4);
% N_T = E_T - Ei;
% dN_T = diff(N_T,T);
% Tmax = vpa(solve(dN_T == 0),4);
% %%Verificação
% % I = -20:0.01:100;
% % N_s = double(subs(N_T,T,I));
% % plot(I,N_s)
% N_max = vpa(subs(N_T,T,Tmax),4);
% N_maxp = N_max/(Emax-Emin)*100;
% %item c
% Eb = vpa(subs(E,[a1 a2 a3 T2 T],[50.37 3.043e-2 -8.567e-5 0 60]),4);
% sa1 = 0.39; sa3 = 1.94e-6; sT2 = 0.019;
% da1 = 60;
% da3 = vpa(3*(-8.567e-5)*60^2,4);
% dT2 = -50.37;
% sE = vpa(sqrt(da1^2*sa1^2+da3^2*sa3^2+da3^2*sT2^2),4);
% 
% nb = round(4.852e-2*Eb);
% sb1 = 0.5/sqrt(3);
% sn = sqrt(sb1^2+(4.852e-2)^2*sE^2);
% 
% Tm = 0.392*151;
% dn = 0.392;
% sErro = sqrt(dn^2*sn^2);

%%% Q2 %%%
% wn = sqrt(0.1024);
% syms csi
% csi = vpa(solve (2*csi*wn == 0.288),4);
% syms s
% G = vpa(0.256/(s^2+2*csi*wn*s+wn^2),4);
% G = collect(simplify(G),s);
% R = 32/s;
% Y = G*R;
% syms t
% DV = collect(simplify(ilaplace(Y,t)),t);
% DV = vpa(simplify(rewrite(DV, 'exp'), 'Steps', 50),4);
% V = DV + 12.5;
% T = 0:0.01:120;
% Vs = double(subs(V,t,T));
% plot(T,Vs)
% G = vpa(0.1024/(s^2+2*csi*wn*s+wn^2),4);
% syms K Ka Kf
% G = G*K*Ka;
% Gmf = K*G*Ka/(1+K*G*Ka*Kf);
% Gmf = vpa(collect(simplify(Gmf),s),4);
% vpa(subs(Gmf,[K Kf],[2.5 0.4]),4);
% G = tf([])