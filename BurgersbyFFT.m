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

dt = 1e-4; %time dimension scale

NoTt = 100; % steps taken for time discretization 

diSx = xL/(NoX-1); % discretized length

visc = 0.1; %  diffusion coefficeint

kx(1:NoX) = 0.;

dkx = 2*pi/(NoX*diSx);

X =0:diSx:2; 

u(1:NoX) = zeros(1,NoX); %preallocating memory

un(1:NoX) = zeros(1,NoX); %preallocating memory

uFFT(1:NoX) = zeros(1,NoX); %preallocating memory

uFFTx(1:NoX) = zeros(1,NoX); %preallocating memory

viscFFT(1:NoX) = zeros(1,NoX); %preallocating memory

viscFFTx(1:NoX) = zeros(1,NoX); %preallocating memory

gauss(1:NoX) = zeros(1,NoX);%preallocating memory

dgauss(1:NoX) = zeros(1,NoX);%preallocating memory




for i=1:NoX/2.
    
         kx(i) = i * dkx; 
end
for i=NoX/2:NoX
    
         kx(i) = (NoX - i) * dkx;
end



for i = 1:NoX    
 
    gauss(i) = exp(-0.25*1i*((X(i)-0)^2/visc));
    dgauss(i) =(-0.5*1i*((X(i)-0)/visc))*exp(-0.25*1i*((X(i)-0)^2/visc));
end

for i = 1:NoX     
    
 u(i) = (-2*visc*(dgauss(i)/gauss(i)))+ 4 ;  % initail Guassian function
 
 %u(i) = (-2*visc*(dgauss(i)))+4;
end


for  it=0:NoTt
    
    un=u;
    h=plot(X,real(u));
    axis([0 2 4 6])
    title({['1-D Burgers'' equation (\nu = ',num2str(visc),')'];['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Transport property profile (u) \rightarrow')
    drawnow;
    refreshdata(h)
           
    
     % fft of diffusion term
         for j = 1:NoX 
        
            viscFFT(j)  = fft(u(j));
    
           %uFFTx = uFFT.*exp(-visc*kx.*kx*dt)+ uFFT.*exp(-1i*kx.*u*dt);  
  
            viscFFTx(j) = viscFFT(j)*(exp(-visc*kx(j)*kx(j)*dt));  
  
            u(j)=ifft(viscFFTx(j));        

         end
         
         %fft of wave stepping term
         
         for j = 1:NoX 
             
           uFFT(j)  = fft(u(j));
    
           %uFFTx = uFFT.*exp(-visc*kx.*kx*dt)+ uFFT.*exp(-1i*kx.*u*dt);  
  
            uFFTx(j) = uFFT(j)*(exp(-1i*kx(j)*u(j)*dt));  
  
            u(j)=ifft(uFFTx(j)); 
            
         end
         
        
 end

