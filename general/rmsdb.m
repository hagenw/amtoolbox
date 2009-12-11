function y = rmsdb(insig)
%RMSDB RMS value of signal (in dB)
%   Usage: y = rmsdb(insig);
%
%   RMSDB(insig) computes the RMS (Root Mean Square) value in dB of a finite 
%   sampled signal sampled at a uniform sampling rate.
%
%   The RMS value of a signal insig of length N is computed by
%
%C                        1  N
%C     rms(insig) = sqrt( - sum insig(n)^2 )
%C                        N n=1
%
%   The level in dB (100 dB for rms(insig) = 1) is calculated by
%
%C 	level(insig) = 20*log10( sqrt( rms(insig).^2 - mean(insig).^2 ) * 10^5 )
%
%   Note: this level is only for the AC component of your signal and the 
%   scaling factor is due to 20*log10(10^5) = 100*log10(10) = 100.
%
%   See also: rms

%   AUTHOR : Hagen Wierstorf
  
% ------ Checking of input parameters ---------

error(nargchk(1,1,nargin));

if ~isnumeric(insig) || ~isvector(insig)
  error('RMS: Input must be a vector.');
end;

% ------ Computation --------------------------

% With var(insig) = rms(insig).^2 - mean(insig).^2 we can 
% compute the RMS value in dB with 
% the following formula:
% y = 20*log10( sqrt(var(insig))*10^5 )
y = 10*log10( var(insig)*10^10 );
