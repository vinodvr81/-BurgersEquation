clear all
clc

tic

dim = 10000;

% a=0.33; %Rabbit growth rate
% b=0.011; %Rabbit death coefficient
% c=0.008; %Fox growth coefficient
% d=0.85; %Fox death rate
% 
% 
% PopN=10000;
% 
% pt=0.4;

a=0.4; %Rabbit growth rate
b=0.001; %Rabbit death coefficient
c=0.001; %Fox growth coefficient
d=0.9; %Fox death rate
% Variables used as input in the ODE-solver
PopN=1000;
pt=0.4;

% initial value problem

RabbitInial=(1-pt)*PopN;
FoxesInial=pt*PopN;

T  = 0:dim;

x = RabbitInial;

y = FoxesInial;

dy = 0.;

dx = 0.;
dt = 1e-2;
% using simple eulers method

for i = 1:dim

dx(i) = dt*(a*x(i) - b*y(i)*x(i));
dy(i)= dt*(c*x(i)*y(i) - d*y(i));


x(i+1) = x(i)+dx(i);

y(i+1) = y(i)+dy(i);


end


figure(1)

plot(T,x,'g')

hold on 
plot(T,y,'--b')

legend('Rabbits','Foxes')
title('Predator-Prey model ODE')
toc

