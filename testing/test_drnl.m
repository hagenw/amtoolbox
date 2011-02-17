%FIXME: There is a huge deviation between drnl and ref_drnl_1 if the
%sampling frequency is high. This must be due to the way the filters are
%calculted. Try with gspi.

if 1
  siglen=10000;
  fs=30000;  
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

outsig=drnl(insig,fs,'jepsen2008');
outsig_ref_1=ref_drnl_1(insig,fs,'jepsen2008');


% There is a deviation. This deviation happends because we convolve the
% filter coefficients, instead of cascading the filters.
res=outsig-outsig_ref_1;
disp('DRNL vs. REF_DRNL_1')
fprintf('Relative l^2-norm: %f\n',norm(res(:))/norm(outsig));
fprintf('Peak SNR: %f\n',20*log10(norm(res(:),inf)/norm(outsig,inf)));

disp('REF_DRNL_1 vs. REF_DRNL')
res=outsig_ref-outsig_ref_1;
fprintf('Relative l^2-norm: %f\n',norm(res(:))/norm(outsig_ref));
fprintf('Peak SNR: %f\n',20*log10(norm(res(:),inf)/norm(outsig_ref,inf)));
