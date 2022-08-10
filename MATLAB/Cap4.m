clear
close all
clc


%%%%%%%%% Ex. 4.1 %%%%%%%%%
% M = 5e-2;
% C = 0.2;
% Ct = M*C;      %Capacitância Térmica
% U1 = 1;
% A = 1e-3;
% Rt1 = 1/(U1*A); %Resistência Térmica
% T1 = Ct*Rt1;
% U2 = 0.2;
% Rt2 = 1/(U2*A); %Resistência Térmica
% T2 = Ct*Rt2;
% 
% t = 0:0.01:300;
% h1 = heaviside(t);
% h2 = heaviside(t-60);
% Ti = 20;
% Tt = 20+80*(1-exp(-60/T1));
% y1 = 20 + 80*(1-exp(-t/T1)).*(h1-h2) -20*h2;
% y2 = (Tt+(Ti-Tt)*(1-exp(-(t-60)/T2))).*h2;
% y = y1 + y2;
% %plot(t,y);
% 
% T = [10 20 50 120 300];
% yv = [100 100 100 20 20];
% for i = 1:5
%     y1 = 20 + 80*(1-exp(-T(i)/T1))*(1-heaviside(T(i)-60)) -20*heaviside(T(i)-60);
%     y2 = (Tt+(Ti-Tt)*(1-exp(-(T(i)-60)/T2)))*heaviside(T(i)-60);
%     y = y1 + y2;
%     fprintf('Erro no tempo %.2f: %.2f\n',T(i),y-yv(i))
% end


%%%%%%%%% Ex. 4.2 %%%%%%%%%

m = 0.5; k = 2e2; l = 6; 
[csi,wn] = param_mola(k,m,l)
syms s
G = 1/k/(1/wn^2*s^2 + 2*csi/wn*s + 1);
G = collect(simplify(G),s);
syms t
collect(simplify(ilaplace(1/s*G,t)),t)

%%%%%%%%% Ex. 4.3 %%%%%%%%%

% csi = 0.1;
% wn = 40;
% syms t s
% w = [10 30 50];
% Fin = 50*[sin(w(1)*t) 10/w(2)*sin(w(2)*t) 10/w(3)*sin(w(3)*t)];
% e = 0;
% for n=1:3
%     [~,pha,mod] = segunda_ordem(csi,wn,w(n));
%     Fc = mod*sin(w(n)*t+rad2deg(pha));
%     e = e + Fc - Fin(n);
% end    
% vpa(e,4);


%%%%%%%%% Ex. 4.4 %%%%%%%%%

% T = 10;
% syms s
% s = tf('s');
% G = 1/(1+T*s);
% 
% % Item a)
% w = largura_banda(G); %%% w = 1/T (Cálc desnec.)
% syms s
% norm(subs(1/(1+T*s),s,w*1i)) ;
% % Item b)
% w = sqrt(((0.95^-1)^2-1)/T^2);
% % Item c)
% T1 = T; T2 = 1;eps=.05
% Gc = tf([T1 1],[T2 1])
% w = freq_max(G*Gc,1-eps)


%%%%%%%%% Ex. 4.5 %%%%%%%%%
% syms s
% s = tf('s');
% k = 10; wn = 10; csi =7; Kd = 1; Ka = 40; R = 350; Kf = 25;
% G1 = 1/k/(1/wn^2*s^2 + 2*csi/wn*s + 1)
% G2 = Kd*Ka*R;
% H = Kf;
% G_mf  = feedback(G1*G2,H)


%%%%%%%%% Ex. 4.6 %%%%%%%%%

% k = 1e2;     %Cte da Mola   
% m = 1;       %Massa do conjunto massa + massa externa
% lambda = 2;  %Fator de Amortecimento
% Krp = 10;    %Sensibilidade - Del_O/Del_T
% g = 9.81;
% m_ext = .5;
% % Item a)
% [csi,wn] = param_mola(k,m,lambda);
% syms s t
% G = Krp/k/(1/wn^2*s^2 + 2*csi/wn*s + 1);
% Fin = m_ext*g/s;
% exp = ilaplace(G*Fin,s,t);
% vpa(exp,2)
% % Item b)
% T = 0:0.1:10;
% y = subs(exp,t,T);
% plot(T,y)

%%%%%%%%% Ex. 4.7 %%%%%%%%%
% wn = 50; csi = 0.2; T = 0.1;
% K = 20;
% Amp = 1e-3*tf([T 0],[T 1]);
% G = tf(50,[1/wn^2 2*csi/wn 1]);
% Gma = zpk(K*Amp*G);
% syms t s
% w = [10 30 50];
% Fin = 50*[sin(w(1)*t) 10/w(2)*sin(w(2)*t) 10/w(3)*sin(w(3)*t)];
% e = 0;
% for n=1:3
%     mod = norm(evalfr(Gma,w(n)*1i));
%     pha = phase(evalfr(Gma,w(n)*1i));
%     Fc = mod*sin(w(n)*t+rad2deg(pha));
%     e = e + Fc - Fin(n);
% end    
% vpa(e,4);


%%%%%%%%% Ex. 4.8 %%%%%%%%%

% M = 5e-2;
% C = 0.2;
% Ct = M*C;      %Capacitância Térmica
% U1 = 1;
% A = 1e-3;
% Rt1 = 1/(U1*A); %Resistência Térmica
% T1 = Ct*Rt1;
% U2 = 0.2;
% Rt2 = 1/(U2*A); %Resistência Térmica
% T2 = Ct*Rt2;
% 
% t = 0:0.01:300;
% h1 = heaviside(t);
% h2 = heaviside(t-60);
% Ti = 50;
% Tt = 50 + 100*(1-exp(-60/T1));
% y1 = 50 + 100*(1-exp(-t/T1)).*(h1-h2) -50.*h2;
% y2 = (Tt+(Ti-Tt)*(1-exp(-(t-60)/T2))).*h2;
% y = y1 + y2;
% plot(t,y)


%%%%%%%%% Ex. 4.9 %%%%%%%%%
% syms t
% T = 63; w0 = 2*pi/T; wc = 2*pi*0.05;
% G = tf(1,[5 1]);
% w = [1 2 3];
% Fin = 10*(sin(w(1)*w0*t) + 1/w(2)*sin(w(2)*w0*t) + 1/w(3)*sin(w(3)*w0*t));
% Fc = 0;
% for n=1:3
%     mod = norm(evalfr(G,w(n)*w0*1i));
%     pha = phase(evalfr(G,w(n)*w0*1i));
%     Fc = Fc + mod*sin(w(n)*w0*t+rad2deg(pha));
% end
% vpa(Fc,4)
% vpa(Fc - Fin,4)


function  wq = freq_max(Gw,Amp)
    [mag,~,wout] = bode(Gw);
    mag = squeeze(mag);
    wq = interp1(20*log10(mag),wout,20*log10(Amp));
end

function  wq = largura_banda(Gw)
    [mag,~,wout] = bode(Gw);
    mag = squeeze(mag);
    wq = interp1(20*log10(mag),wout,20*log10(2^-0.5));
end

function [G,pha,mod] = segunda_ordem(csi,wn,w)
    syms s
    G = 1/(s^2/wn^2+2*s*csi/wn+1);
    pha = phase(subs(G,s,w*1i));
    mod = norm(subs(G,s,w*1i));
end

function [csi,wn] = param_mola(k,m,lambda)
    wn = sqrt(k/m);
    csi = lambda/(2*sqrt(k*m));
end
    