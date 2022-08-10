clear
close all
clc 

%%% Q1 %%%
syms x
eqn = 1/(1+100/140*x) - 1/(1+x) == 0.1/15;
vpa(solve(eqn),4)