function outsig=filterbankz(b,a,insig,hopsize)
%FILTERBANKZ  Filter bank with zero boundary condition
%   Usage: outsig=filterbankz(b,a,insig);
%          outsig=filterbankz(b,a,insig,hopsize);
%
%   `filterbankz(b,a,insig,hopsize)` filters the input signal with the
%   filters described in *a* and *b*. *hopsize* is a vector with a length
%   equal to the number of filters. Each channel is sub-sampled by the
%   corresponding hopsize.
%
%   If *a* and *b* are matrices then each row corresponds to a subband
%   channel.
%
%   If *insig* is a matrix then filtering is applied along the columns.
%
%   The output coefficients are stored a cell array. More precisely, the
%   n'th cell of *c*, `c{m}`, is a 2D matrix of size $M(n) \times W$ and
%   containing the output from the m'th channel subsampled at a rate of
%   $a(m)$.  `c{m}(n,l)` is thus the value of the coefficient for time index
%   *n*, frequency index *m* and signal channel *l*.
%
%   See also: gammatone, ufilterbankz, auditoryfilterbank
  
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

