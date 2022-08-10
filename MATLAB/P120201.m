clear
close all
clc


%%% Q1 %%%
% syms T
% Eideal = 5/120*T;
% N = 6.933e-2*T - 2.305e-4*T^2 - Eideal;
% dN = diff(N,T);
% TNmax = vpa(solve(dN == 0),4);
% vpa(subs(N,T,TNmax),4)/5;
% 
% syms C0 C1 C2
% E = C0 + C1*T + C2*T^2
% dC0 = 3.21e-3;dC1=4.83e-3;
% delC0 = 1
% delC1 = 6.933e-2
% delC2 = vpa(2*(- 2.305e-4)*TNmax,4)
% dE = sqrt(dC0^2*delC0^2+dC1^2*delC1^2)


%%% Q2 %%%

% syms s K Kf Ka T
% G = 1/(T*s+1);
% Gmf = K*G*Ka/(1+K*G*Ka*Kf);
% Gmf = collect(simplify(Gmf),s);
% 
% wB1 = (K*Ka*Kf+1)/T;
% wB1 = subs(wB1,[T Kf],[1e-3 1/Ka]);
% vpa(solve(wB1 == 1e4),4)