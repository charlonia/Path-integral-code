clc
clear
grid_div=0.004;
time_div=0.002;
m=1;
h=1;
z = i/(2*time_div)*m/h
f = @(x) exp(z.*(x.^2));
x=1;
A=0:grid_div:x;
G = f(A)*grid_div;
plot(A,real(G));
g=2*sum(G)*sqrt(m/(2*pi*i*h*time_div))
