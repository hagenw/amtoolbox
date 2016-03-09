function [outsig,fb]=filterbank_block(insig,fb)
%FILTERBANK_BLOCK  Process signal block in a filterbank
%   Usage: [outsig,fb]=filterbank(fb,insig);
%
%   `[outsig,fb]=filterbank_block(insig,fb)` filters a block of input using the
%   filters described in the structrure *fb*. The structure *fb* must first
%   have been created using |filterbank_init|.
%
%   If *insig* is a matrix then filtering is applied along the columns. In
%   this case, the output will be 3-dimensional. First dimension
%   is time, second dimension is frequency channel and third dimension
%   corresponds to the columns of the input signal.
%
%   See also: filterbank_init, ufilterbankz
  
%   AUTHOR : Peter L. Søndergaard

% ------ Checking of input parameters ---------  

error(nargchk(2,2,nargin));


% ------ Computation --------------------------

siglen=size(insig,1);

if rem(siglen,fb.hopsize)~=0
  error(['FILTERBANK_BLOCK: The hopsize must divide the length of the ' ...
         'input block.']);
end;

fb.outlen=ceil(siglen/fb.hopsize);

outsig=zeros(fb.outlen,fb.nchannels,fb.nsigs);

for ii=1:fb.nchannels
  [res,fb.zi(:,:,ii)] = filter(fb.b(ii,:),fb.a(ii,:),insig,squeeze(fb.zi(:,:,ii)));
  res = res(1:fb.hopsize:siglen,:);  
  outsig(:,ii,:) = reshape(res,fb.outlen,1,fb.nsigs);
end;

fb.outstart = fb.outend+1;
fb.outend   = fb.outstart+fb.outlen-1;

