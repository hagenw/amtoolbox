function y = rms(x)
%RMS RMS value of signal
%   Usage: y = rms(x)
%
%   RMS(x) computes the RMS (Root Mean Square) value of a finte sampled
%   signal sampled at a uniform sampling rate.
%
%   The RMS value of a signal x of length N is computed by
%M
%C                    1  M
%C     rms(x) = sqrt( - sum x(n)^2 )
%C                    N n=1
%M

%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ---------
  
error(nargchk(1,1,nargin));

if ~isnumeric(x) || ~isvector(x)
  error('RMS: Input must be a vector.');
end;

% ------ Computation --------------------------

% It is better to use 'norm' instead of explicitly summing the squares, as
% norm (hopefully) attempts to avoid numerical overflow. Matlab does not
% have 'sumsq'.
y = norm(x)/sqrt(length(x));
