t=0:0.01:10;
y=sin(t);
%plot(t,y);
%-------------------------
size(y)
size(t)
dy=diff(y)./diff(t)