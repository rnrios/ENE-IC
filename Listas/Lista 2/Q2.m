clear
close all
clc


syms k b

eq1 = [100 == k*exp(b/218), 1e4 == k*exp(b/393)];
x = solve(eq1, 'Real', true);
b = vpa(x.b,4);
k = vpa(x.k,4);

Rmid = k*exp(b/305.5)
Rmax = 1e4;
Rmin = 1e2;

syms mid min max r4
eq = 2*mid/(mid+r4) == min/(min+r4)+max/(max+r4);
r4 = solve(eq,r4);
r4 = r4(2)

r4 = subs(r4, [max, mid, min], [Rmax, Rmid, Rmin])
r = r4/Rmin
Vs = 1/(Rmax/(r4+Rmax)-1/(1+r))