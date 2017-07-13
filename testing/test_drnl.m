function test_failed=test_drnl
%TEST_DRNL  Test the DRNL
%
%  Test if the DRNL produce the same results as the reference implementation.
%
%  

amt_disp(' ===============  TEST_DRNL ================');

% The reference implementations are not dbspl-safe, so switch to the
% dboffset=100 setting.
  
save_dboffset = dbspl(1);

ltfatsetdefaults('dbspl','dboffset',100);

testsig=2;

switch testsig
  case 1
   siglen=10000;
   fs=16000;  
   insig=randn(siglen,1);
 case 2
  insig=greasy;
  fs=16000;
  siglen=length(insig);    
 case 3
  insig=gspi;
  fs=44100;
  siglen=length(insig);
 case 4
  siglen=20000;
  fs=441000;  
  insig=randn(siglen,1);
end;

fc=erbspacebw(80,8000);
nfc=length(fc);

outsig_ref=zeros(siglen,nfc);

for ii=1:nfc
  
  outsig_ref(:,ii)=ref_drnl(insig,fc(ii),fs);
end;

outsig=drnl(insig,fs,'jepsen2008');
outsig=gaindb(outsig,50);

outsig_ref_1=ref_drnl_1(insig,fs,'jepsen2008');


% There is a deviation. This deviation happends because we convolve the
% filter coefficients, instead of cascading the filters.
res=outsig-outsig_ref_1;
amt_disp('DRNL vs. REF_DRNL_1')
fprintf('Relative l^2-norm: %f\n',norm(res(:))/norm(outsig(:)));
fprintf('Peak SNR: %f\n',20*log10(norm(res(:),inf)/norm(outsig(:),inf)));

amt_disp('REF_DRNL_1 vs. REF_DRNL')
res=outsig_ref-outsig_ref_1;
fprintf('Relative l^2-norm: %f\n',norm(res(:))/norm(outsig_ref(:)));
fprintf('Peak SNR: %f\n',20*log10(norm(res(:),inf)/norm(outsig_ref(:),inf)));

% Test binaural
binoutsig=drnl([insig,insig],fs,'jepsen2008');
binoutsig=gaindb(binoutsig,50);

res=2*outsig-binoutsig(:,:,1)-binoutsig(:,:,2);
fprintf('Binaural: %f\n',norm(res(:)));

ltfatsetdefaults('dbspl','dboffset',save_dboffset);


