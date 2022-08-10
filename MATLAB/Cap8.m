clear;
close all;
clc;

%%% Q1 %%%
% syms a b
% R0 = 100;
% X = solve(R0*(1+a*100+b*100^2) == 138.5, R0*(1+a*200+b*200^2) == 175.83)
% vpa(X.a,4)
% vpa(X.b,4)
 
%%% Q2 %%%
l=0.25; w=.06; t=3e-3; E=7e10; G=2.1; R0=120; F=0.5;
x = l/2;
e = 6*(l-x)*F/(w*t^2*E);
R1 = R0*(1+G*e)
R2 = R0*(1-G*e)
