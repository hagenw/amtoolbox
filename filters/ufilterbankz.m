function outsig=ufilterbankz(b,a,insig,hopsize)
%UFILTERBANKZ  Uniform Filter bank with zero boundary condition
%   Usage: outsig=ufilterbankz(b,a,insig);
%          outsig=ufilterbankz(b,a,insig,hopsize);
%
%   `ufilterbankz(b,a,insig)` filters the input signal with the filters
%   described in *a* and *b*.
%
%   `ufilterbankz(b,a,insig,hopsize)` does the same, but only outputs every
%   *hopsize* sample in the time domain.
%
%   If *a* and *b* are matrices then each row corresponds to a subband
%   channel.
%
%   If *insig* is a matrix then filtering is applied along the columns.
%
%   If both *insig*, *a* and *b* are matrices the output will be
%   3-dimensional. First dimension is time, second dimension is frequency
%   and the third dimension corresponds to the columns of the input signal.
%
%   See also: gammatone, filterbankz, auditoryfilterbank
  
%   AUTHOR : Peter L. Soendergaard

%% ------ Checking of input parameters ---------  

error(nargchk(3,4,nargin));

if nargin==3
  hopsize=1;
end;


%% ------ Computation --------------------------

[insig,siglen,dummy,nsigs,dim,permutedsize,order]=assert_sigreshape_pre(insig,[],[], ...
                                                  upper(mfilename));
nchannels=size(b,1);

outlen=ceil(siglen/hopsize);

outsig=zeros(outlen,nchannels,nsigs);

for ii=1:nchannels
  res = filter(b(ii,:),a(ii,:),insig);
  res = res(1:hopsize:siglen,:);  
  outsig(:,ii,:) = reshape(res,outlen,1,nsigs);
end;

% Modify permutedsize and order to reflect that first dimension has split
% in two
permutedsize=[outlen,nchannels,permutedsize(2:end)];
order=assert_groworder(order);

outsig=assert_sigreshape_post(outsig,dim,permutedsize,order);

