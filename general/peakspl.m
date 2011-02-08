function y = peakspl(insig,options)
%PEAKSPL Peak value of signal (in dB)
%   Usage: y = peakspl(insig);
%          y = peakspl(insig,'ac');
%
%   PEAKSPL(insig) computes the peak SPL value in dB SPL of a finite sampled
%   signal sampled at a uniform sampling rate. The dB value is computed from
%   the convention that a pure tone at 100 dB SPL has an RMS value of 1.
%
%   PEAKSPL(x,'ac') does the same, but considers only the AC component of the
%   signal (i.e. the mean is removed).
%
%   See also: gaindb, crestfactor, setleveldb, rmsdb
%
%R  moore2003introduction

%   AUTHOR : Hagen Wierstorf
  
% ------ Checking of input parameters ---------

error(nargchk(1,2,nargin));

if ~isnumeric(insig) || ~isvector(insig)
  error('PEAKSPL: Input must be a vector.');
end;

if (nargin==1) || (~ischar(options))
  options='';
end;

% ------ Computation --------------------------

% XXX Is this correct???
y = 20*log10( norm(insig,Inf))+100;
