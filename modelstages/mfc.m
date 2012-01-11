function mfc=mfc(fc)

Q = 2;
bw = 5;
ex=(1+1/(2*Q))/(1-1/(2*Q));

startmf = 5;

umf = min(fc.*0.25, 1000);  

tmp = fix((min(umf,10) - startmf)/bw);
tmp = 0:tmp;
mfc = startmf + 5*tmp;
tmp2 = (mfc(end)+bw/2)/(1-1/(2*Q));
tmp = fix(log(umf/tmp2)/log(ex));
tmp = 0:tmp;
tmp = ex.^tmp;
mfc=[0 mfc tmp2*tmp];

%OLDFORMAT
