function f=comp_ifftreal(c,N);
%COMP_IFFTREAL  Inverse FFT for real valued signals.
%   Usage: f=comp_ifftreal(c,N);
%

%   AUTHOR : Peter Soendergaard

W=size(c,2);
N2=size(c,1);
  
cfull=zeros(N,W);
cfull(1:N2,:)=c;

if rem(N,2)==0
  cfull(N2+1:N,:)=conj(c(N2-1:-1:2));
else
  cfull(N2+1:N,:)=conj(c(N2:-1:2));
end;

% Force IFFT along dimension 1, since we have permuted the dimensions
% manually
f=ifft(cfull,N,1);
