if 0
  siglen=10000;
  fs=16000;  
  insig=randn(siglen,1);

else
  insig=greasy;
  fs=16000;
  siglen=length(insig);    
end;

fc=erbspacebw(80,8000);
nfc=length(fc);

outsig_ref=zeros(siglen,nfc);

for ii=1:nfc
  
  outsig_ref(:,ii)=ref_drnl(insig,fc(ii),fs);
end;

outsig=drnl(insig,fc,fs);

res=outsig-outsig_ref;

% There is a deviation. This deviation happends because we convolve the
% filter coefficients, instead of cascading the filters.
norm(res(:))/norm(outsig)