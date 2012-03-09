function fb=filterbank_init(b,a,varargin);
%FILTERBANK_INIT  Wrapper around filter to multiple filters
%   Usage: outsig=filterbank_init(b,a);
%          outsig=filterbank_init(b,a,nsigs,hopsize);
%
%   `fb=filterbank_init(b,a)` creates a filterbank structure *fb* for use with
%   |filterbank_block|_. The filterbank will filter the input signals with the
%   filters described in *b* and *a*.
%
%   `fb=filterbank_init(b,a,nsigs)` does the same assuming that the input
%   to |filterbank_block|_ will consist of *nsigs* signal at once.
%
%   See also: ufilterbankz, filterbank_block

%   AUTHOR : Peter L. Soendergaard

% ------ Checking of input parameters ---------  

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

definput.keyvals.nsigs=1;
definput.keyvals.hopsize=1;

[flags,keyvals,fb.nsigs,fb.hopsize]  = ltfatarghelper({'nsigs','hopsize'},definput,varargin);

zilen=max(size(a,2),size(b,2))-1;

fb.b=b;
fb.a=a;
fb.nchannels=size(b,1);

fb.outstart = 0;
fb.outlen   = 0;
fb.outend   = 0;

% Initialize the initial conditions to zero.
fb.zi=zeros(zilen,fb.nsigs,fb.nchannels);
