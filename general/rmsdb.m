function y = rmsdb(insig)
%RMSDB RMS value of signal (in dB)
%   Usage: y = rmsdb(insig);
%
%   RMSDB(insig) computes the RMS (Root Mean Square) value in dB SPL of a 
%   finite sampled signal sampled at a uniform sampling rate. The dB value is
%   computed from the convention that 100 dB SPL corresponds to an RMS
%   value of 1.
%
%   See also: rms, gaindb, setleveldb

%   AUTHOR : Hagen Wierstorf
  
% ------ Checking of input parameters ---------

error(nargchk(1,1,nargin));

if ~isnumeric(insig) || ~isvector(insig)
  error('RMSDB: Input must be a vector.');
end;

% ------ Computation --------------------------

% The level of a signal in dB SPL is given by the following formula:
% level = 20*log10(p/p_0)
y = 20*log10( rms(insig) )+100;
