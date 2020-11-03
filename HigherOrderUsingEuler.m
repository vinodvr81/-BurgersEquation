clear all
clc
tic
dim =10000;

t=1:dim;

x = 7;
y =13;

dx = 0.;
dy =0.;

dt = 1e-3;

loop = 1;

while loop < dim

dx(loop) = dt*y(loop);    
dy(loop) = dt*(1/2*(11*exp(-t(loop)^3)-3*x(loop) - 5*x(loop)));


x(loop+1) =  x(loop)+ dx(loop);
y(loop+1) = y(loop)+ dy(loop);
loop = loop+1;

end

figure(1)
plot(x,'g')
hold on 
plot(y,'--b')

toc

