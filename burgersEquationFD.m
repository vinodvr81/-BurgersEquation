
clear
clc
clf

% The Burgers equation can be written as  u_t + u*u_x = v*u_xx

% Wave stepping term is u*u_x

% Diffusion term is v*u_xx

% assuming v = 0

% discretization in space 1m total length with 50 points

% discretization in time 1e-2 sec total length with 50 points

% initialization of parameters

xL = 2; % space dimension scale

NoX = 20; % steps taken for space discretization

dt = 0.01; %time dimension scale

NoTt = 50; % steps taken for time discretization 

diSx = xL/(NoX-1); % discretized length

visc = 0.1; %  diffusion coefficeint

X =0:diSx:2; 

u(1:NoX) = zeros(1,NoX); %preallocating memory
un(1:NoX) = zeros(1,NoX);

ip(1:NoX) = zeros(1,NoX); %preallocating memory
im(1:NoX) = zeros(1,NoX);

gauss(1:NoX) = zeros(1,NoX);%preallocating memory
dgauss(1:NoX) = zeros(1,NoX);%preallocating memory

for i = 1:NoX    
    ip(i)=i+1;
    im(i)=i-1;
    gauss(i) = exp(-0.25*((X(i)-0)^2/visc));
    dgauss(i) =(-0.5*((X(i)-0)/visc))*exp(-0.25*((X(i)-0)^2/visc));
end

ip(NoX)=1;
im(1)=NoX;

%plot(gauss)

for i = 1:NoX     
    
 u(i) = (-2*visc*(dgauss(i)/gauss(i)))+ 4;  % initail Guassian function
 
 %u(i) = (-2*visc*(dgauss(i)))+4;
end


 %plot(u)
for  it=0:NoTt
    un=u;
    h=plot(X,u);
    axis([0 2 4 6])
    title({['1-D Burgers'' equation (\nu = ',num2str(visc),')'];['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Transport property profile (u) \rightarrow')
    drawnow;
    refreshdata(h)
    for j = 1:NoX        
        u(j) =un(j)+ (visc*dt*((un(ip(j)) - 2*un(j)+ un(im(j)))/diSx^2) - (dt*un(j)*((un(j) - un(im(j)))/diSx)));           
    

    end
end












