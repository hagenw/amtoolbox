function [b,a]=gammatone(fc,fs,varargin);
%GAMMATONE  Gammatone filter coefficients
%   Usage: [b,a] = gammatone(fc,fs,n,betamul);
%          [b,a] = gammatone(fc,fs,n);
%          [b,a] = gammatone(fc,fs);
%
%   Input parameters:
%      fc    -  center frequency in Hz.
%      fs    -  sampling rate in Hz.
%      n     -  filter order.
%      beta  -  bandwidth of the filter.
%
%   Output parameters:
%      b     -  nominator coefficients.
%      a     -  denominator coefficients.
%
%   GAMMATONE(fc,fs,n,betamul) computes the filter coefficients of a digital
%   gammatone filter with center frequency fc, order n, sampling rate fs and
%   bandwith determined by betamul. The bandwidth _beta of each filter is
%   determined as betamul times AUDFILTBW of the center frequency of
%   corresponding filter.
%
%   GAMMATONE(fc,fs,n) will do the same but choose a filter bandwidth
%   according to Glasberg and Moore (1990). The order n can only be 2 or
%   4. If the order is 4 then betamul is choosen to be 1.0183 and if the
%   order is 2 then betamul is choosen to be 0.637.
%
%   GAMMATONE(fc,fs) will do as above for a 4th order filter.
%
%   If fc is a vector, each entry of fc is considered as one center
%   frequency, and the corresponding coefficients are returned as row
%   vectors in the output.
%
%   The inpulse response of the gammatone filter is given by
%
%M    g(t) = a*t^(n-1)*cos(2*pi*fc*t)*exp(-2*pi*beta*t)
%F  \[g(t) = at^{n-1}cos(2\pi\cdot fc\cdot t)e^{-2\pi \beta \cdot t}\]
%
%   GAMMATONE( ..., 'complex') will generate filter coefficients
%   corresponding to a complex valued filterbank modulated by exponential
%   functions. This is usefull for envelope extration purposes.
%
%   To create the filter coefficients of a 1-erb spaced filter bank using
%   gammatone filters use the following construction
%
%C    [b,a] = gammatone(erbspacebw(flow,fhigh),fs);
%  
%R  aertsen1980strI glasberg1990daf
  
%   AUTHOR : Stephan Ewert, Peter L. Soendergaard

% ------ Checking of input parameters ---------
  

if nargin<2
  error('Too few input arguments.');
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(fc) || ~isvector(fc) || any(fc<0) || any(fc>fs/2)
  error(['%s: fc must be a vector of positive values that are less than half ' ...
         'the sampling rate.'],upper(mfilename));
end;

defnopos.keyvals.n=4;
defnopos.keyvals.betamul=[];
defnopos.flags.real={'real','complex'};

[flags,keyvals,n,betamul]  = ltfatarghelper({'n','betamul'},defnopos,varargin);

if ~isnumeric(n) || ~isscalar(n) || n<=0 || fix(n)~=n
  error('%s: n must be a positive, integer scalar.',upper(mfilename));
end;

if isempty(betamul)
  switch(n)
   case 2
    betamul =  0.637;
   case 4
    betamul = 1.0183;
   otherwise
    error(['GAMMATONE: Default value for beta can only be computed for 2nd ' ...
           'and 4th order filters.']);
  end;
else
  if ~isnumeric(betamul) || ~isscalar(betamul) || betamul<=0
    error('%s: beta must be a positive scalar.',upper(mfilename));
  end;
end;


% ------ Computation --------------------------

nchannels = length(fc);

b=zeros(nchannels,1);
a=zeros(nchannels,n+1);

% ourbeta is used in order not to mask the beta function.

ourbeta = betamul*audfiltbw(fc);
  
for ii = 1:nchannels
  
  btmp=1-exp(-2*pi*ourbeta(ii)/fs);
  atmp=[1, -exp(-(2*pi*ourbeta(ii) + i*2*pi*fc(ii))/fs)];
  
  b2=1;
  a2=1;
  
  for jj=1:n
    b2=conv(btmp,b2);
    a2=conv(atmp,a2);
  end
  
  % Place the result (a row vector) in the output matrices.
  b(ii,:)=b2;
  a(ii,:)=a2;

end;

if flags.do_real
  b=real(b);
  a=real(a);
end;

