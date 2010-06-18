function y = rmsdb(insig,options)
%RMSDB RMS value of signal (in dB)
%   Usage: y = rmsdb(insig);
%          y = rmsdb(insig,'ac');
%
%   RMSDB(insig) computes the RMS (Root Mean Square) value in dB SPL of a
%   finite sampled signal sampled at a uniform sampling rate. The dB value
%   is computed from the convention that a pure tone at 100 dB SPL has an
%   RMS value of 1.
%
%   RMSDB(x,'ac') does the same, but considers only the AC component of the
%   signal (i.e. the mean is removed).
%
%   See also: rms, gaindb, crestfactor, setleveldb
%
%R  moore2003introduction

%   AUTHOR : Hagen Wierstorf
  
% ------ Checking of input parameters ---------

error(nargchk(1,2,nargin));

if ~isnumeric(insig) || ~isvector(insig)
  error('RMSDB: Input must be a vector.');
end;

if (nargin==1) || (~ischar(options))
  options='';
end;

% ------ Computation --------------------------

% The level of a signal in dB SPL is given by the following formula:
% level = 20*log10(p/p_0)
% To get to the standard used in the toolbox.
y = 20*log10( rms(insig,options) )+100;
