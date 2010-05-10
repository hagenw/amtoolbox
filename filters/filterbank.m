function outsig=filterbank(b,a,insig,hopsize)
%FILTERBANK  Wrapper around filter to multiple filters
%   Usage: outsig=filterbank(b,a,insig);
%          outsig=filterbank(b,a,insig,hopsize);
%
%   FILTERBANK(b,a,insig) filters the input signal with the filters
%   described in _a and b.
%
%   FILTERBANK(b,a,insig,hopsize) does the same, but only outputs every
%   hopsize sample in the time domain.
%
%   If _a and b are matrices then each row corresponds to a subband
%   channel.
%
%   If insig is a matrix then filtering is applied along the columns.
%
%   If both insig, _a and b are matrices the output will be
%   3-dimensional. First dimension is time, second dimension is frequency
%   channel third dimension corresponds to the column of the input signal.

%   AUTHOR : Peter L. Soendergaard

% ------ Checking of input parameters ---------  

error(nargchk(3,4,nargin));

if nargin==3
  hopsize=1;
end;


% ------ Computation --------------------------

nchannels=size(b,1);

siglen=size(insig,1);
nsigs=size(insig,2);

outlen=ceil(siglen/hopsize);

outsig=zeros(outlen,nchannels,nsigs);

for ii=1:nchannels
  res = filter(b(ii,:),a(ii,:),insig);
  res = res(1:hopsize:siglen,:);  
  outsig(:,ii,:) = reshape(res,outlen,1,nsigs);
end;
