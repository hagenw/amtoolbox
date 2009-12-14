function y = rms(insig,options)
%RMS RMS value of signal
%   Usage: y = rms(insig);
%          y = rms(insig,'ac');
%
%   RMS(x) computes the RMS (Root Mean Square) value of a finite sampled
%   signal sampled at a uniform sampling rate.
%
%   RMS(x,'ac') does the same, but considers only the AC component of the
%   signal (i.e. the mean is removed).
%
%   The RMS value of a signal x of length N is computed by
%
%C                    1  N
%C     rms(x) = sqrt( - sum x(n)^2 )
%C                    N n=1
%
%R  moore2003introduction

%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ---------
  
error(nargchk(1,2,nargin));

if ~isnumeric(insig) || ~isvector(insig)
  error('RMS: Input must be a vector.');
end;

if (nargin==1) || (~ischar(options))
  options='';
end;

% ------ Computation --------------------------

% It is better to use 'norm' instead of explicitly summing the squares, as
% norm (hopefully) attempts to avoid numerical overflow.

switch(lower(options))
 case 'ac'
  y = norm(insig-mean(insig))/sqrt(length(insig));
 otherwise
  y = norm(insig)/sqrt(length(insig));
end;