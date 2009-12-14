function y = rms(insig)
%RMS RMS value of signal
%   Usage: y = rms(insig)
%
%   RMS(insig) computes the RMS (Root Mean Square) value of a finite sampled
%   signal sampled at a uniform sampling rate.
%
%   The RMS value of a signal insig of length N is computed by
%M
%C                    1  N
%C     rms(insig) = sqrt( - sum insig(n)^2 )
%C                    N n=1
%M
%
%R  moore2003introduction

%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ---------
  
error(nargchk(1,1,nargin));

if ~isnumeric(insig) || ~isvector(insig)
  error('RMS: Input must be a vector.');
end;

% ------ Computation --------------------------

% It is better to use 'norm' instead of explicitly summing the squares, as
% norm (hopefully) attempts to avoid numerical overflow. Matlab does not
% have 'sumsq'.
y = norm(insig)/sqrt(length(insig));
