function y = dbspl(insig,options)
%DBSPL RMS value of signal (in dB)
%   Usage: y = dbspl(insig);
%          y = dbspl(insig,'ac');
%
%   DBSPL(insig) computes the SPL (sound pressure level) of the input signal
%   measured in dB, using the convention that a pure tone at 100 dB SPL has
%   an RMS value of 1.
%
%   DBSPL(x,'ac') does the same, but considers only the AC component of the
%   signal (i.e. the mean is removed).
%
%   See also: setleveldb
%
%R  moore2003introduction

%   AUTHOR : Hagen Wierstorf
  
% ------ Checking of input parameters ---------

error(nargchk(1,2,nargin));

if ~isnumeric(insig) || ~isvector(insig)
  error('DBSPL: Input must be a vector.');
end;

if (nargin==1) || (~ischar(options))
  options='';
end;

% ------ Computation --------------------------

% The level of a signal in dB SPL is given by the following formula:
% level = 20*log10(p/p_0)
% To get to the standard used in the toolbox.
y = 20*log10( rms(insig,options) )+100;
