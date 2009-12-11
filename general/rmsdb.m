function y = rmsdb(insig)
%RMSDB RMS value of signal (in dB)
%   Usage: y = rmsdb(insig);
%
%   RMSDB(insig) computes the RMS (Root Mean Square) value in dB SPL of a finite 
%   sampled signal sampled at a uniform sampling rate. The Db value is
%   computed from the convention that 100 Db SPL corresponds to an RMS
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

% With var(insig) = rms(insig).^2 - mean(insig).^2 we can 
% compute the RMS value in dB with 
% the following formula:
% y = 20*log10( sqrt(var(insig))*10^5 )
y = 20*log10( rms(insig))+100;
