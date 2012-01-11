% Usage: [b,a]=getGFBFilterCoef(nrch,centerf,bwerb,fs)

function [b,a]=getGFBFilterCoefs(nrch,centerf,bwerb,fs);

b=[];
a=[];

%centerf=24.7*9.265*(exp(centererb/9.265)-1);
erb = 24.7 + centerf/9.265;
beta=1.0183*erb*bwerb;     % 4th order gamma, 0.637 if 2nd order.

for i=1:nrch
   [bt,at]=gammaf1(centerf(i),beta(i),4,2,fs);
   b=[b;bt];
   a=[a;at];
end

%OLDFORMAT
