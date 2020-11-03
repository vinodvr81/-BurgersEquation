%step 1: clears memory and graphics window 
clear
clc
clf

tic
%step 2:initialization of parameters
psatur =   10;
pb1 = 0.05;
pb2 = 0.005;

diffusioncoeff= 1.5e-4;
dt = 1e-3;
palpha = 2;
pbeta = 0;

trans = 0.1;
propagationlength = 1;

numberofturns = 1; 
dimx = 16; 
dimz = 16;

chipsize = 20;
dx = chipsize/(dimx-1) ;
dz = propagationlength/(dimz-1);
piagr(1:dimx,1:dimz) = 1.44;
piasa(1:dimx,1:dimz) = 0.5;

DcarG(1:dimx,1:dimz) =0.;
dcarA(1:dimx,1:dimz) = 0.;
ef(1:dimx,1:dimz) =(0+0*1i);

 

for j = 1:dimx
for i = 1:(dimz/2-3)    
ef(j,i) = (0.05+0.05*1i);
DcarG(j,i) =  1.15; 
dcarA(j,i) = -0.57;

end
for i = (dimz/2-2):(dimz/2+2)
ef(j,i) = (0.3+0.3*1i);
DcarG(j,i) =  1.43; 
dcarA(j,i) = -0.23;

end
for i = (dimz/2+3):dimz
ef(j,i) = (0.05+0.05*1i);
DcarG(j,i) =  1.15; 
dcarA(j,i) = -0.57;

end
end




noise = 0.0001;


 
 while numberofturns < 600000
 
    def = dt*(((1-1i*palpha).*DcarG+  (1-1i*pbeta).*dcarA-1).*ef + noise*rand(dimx,dimz));         
    dDcarG = -dt*pb1*((DcarG.*(1 + abs(ef).^2))- piagr);       
    ddcarA = -dt*pb2*((dcarA.*(1 + psatur*abs(ef).^2))+ piasa); 
    
    ef = def + ef;
    DcarG = DcarG + dDcarG;
    dcarA = dcarA + ddcarA;   
  
   numberofturns = numberofturns + 1;
   
  power(numberofturns) =  sum(sum(abs((ef).^2)));   
  
   end
     
   
 plot (power)
   
 toc
   
   
   
   
   
   
   
   
   
   
   
   
   
   