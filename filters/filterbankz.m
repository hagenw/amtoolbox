function outsig=filterbankz(b,a,insig,hopsize)
%FILTERBANKZ  Filter bank with zero boundary condition
%   Usage: outsig=filterbank(b,a,insig);
%          outsig=filterbank(b,a,insig,hopsize);
%
%   FILTERBANKZ(b,a,insig,hopsize) filters the input signal with the filters
%   described in _a and b. hopsize is a vector with a length equal to the number
%   of filters. Each channel is sub-sampled by the correcsponding hopsize.
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

%% ------ Checking of input parameters ---------  

error(nargchk(4,4,nargin));


%% ------ Computation --------------------------

[insig,siglen,dummy,nsigs,dim,permutedsize,order]=assert_sigreshape_pre(insig,[],[], ...
                                                  upper(mfilename));
nchannels=size(b,1);


outsig=cell(nchannels,1);

for ii=1:nchannels
  % Calculate the new length in the time domain of this channel
  outlen=ceil(siglen/hopsize);

  % Do the actual filtering.
  res = filter(b(ii,:),a(ii,:),insig);

  % Subsample the output, reshape a multidimensional array to the correct size and store.
  permutedsize(1)=outlen;
  outsig{ii} = assert_sigreshape_post(res(1:hopsize:siglen,:),dim,permutedsize,order);  
end;

