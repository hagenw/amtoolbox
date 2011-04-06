function outsig = noise(siglen,nsigs,varargin)
% NOISE Generates a noise signal
%   Usage: outsig = noise(siglen,nsigs,type);
%
%   Input parameters:
%       siglen    - Length of the noise (samples)
%       nsigs     - Number of signals (default is 1)
%       type      - type of noise: 'white', 'brown', 'pink', 'red'
%
%   Output parameters:
%       outsig      - siglen x nsigs signal vector
%
%   NOISE(siglen,nsigs,type) generates nsigs channels containing noise of the
%   given type with the length of siglen. The signals are arranged as columns
%   in the output.
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameter -------------------------------------
nargmin = 1;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(siglen) || ~isscalar(siglen) || siglen<=0
    error('%s: siglen has to be a positive scalar.',upper(mfilename));
end

if nargin<2
    nsigs = 1;
elseif ~isnumeric(nsigs) || ~isscalar(nsigs) || nsigs<=0
    error('%s: nsigs has to be a positive scalar.',upper(mfilename));
end

defnopos.flags.real={'white','pink','brown','red'};

[flags,keyvals]  = ltfatarghelper({},defnopos,varargin);

if flags.do_white
  outsig=randn(siglen,nsigs);
end;

if flags.do_brown || flags.do_red
  outsig=cumsum(randn(siglen,nsigs));
end;

if flags.do_pink
  % --- Handle trivial condition

  if siglen==1
    outsig=ones(1,nsigs);
    return;
  end;

  % ------ Computation -----------------------------------------------------
  fmax = floor(siglen/2)-1;
  f = (2:(fmax+1)).';
  % 1/f amplitude factor
  a = 1./sqrt(f);
  % Random phase
  p = randn(fmax,nsigs) + i*randn(fmax,nsigs);
  sig = repmat(a,1,nsigs).*p;

  outsig = ifftreal([ones(1,nsigs); sig; 1/(fmax+2)*ones(1,nsigs)],siglen);

end;

% Scale output
outsig = outsig ./ (max(abs(outsig(:)))+eps);

