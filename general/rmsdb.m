function y = rmsdb(x)
%RMSDB RMS value of signal (in dB)
%   Usage: y = rmsdb(x);
%
%   RMSDB(x) computes the RMS (Root Mean Square) value in dB of a finite 
%   sampled signal sampled at a uniform sampling rate.
%
%   The RMS value of a signal x of length N is computed by
%
%C                    1  M
%C     rms(x) = sqrt( - sum x(n)^2 )
%C                    N n=1
%
%   The level in dB (100 dB for rms(x) = 1) is calculated by
%
%C 	level(x) = 20*log10( sqrt( rms(x).^2 - mean(x).^2 ) * 10^5 )
%
%   Note: this level is only for the AC component of your signal and the 
%   scaling factor is due to 20*log10(10^5) = 100*log10(10) = 100.
%
%   See also: rms

%   AUTHOR : Hagen Wierstorf
  
% ------ Checking of input parameters ---------

error(nargchk(1,1,nargin));

if ~isnumeric(x) || ~isvector(x)
  error('RMS: Input must be a vector.');
end;

% ------ Computation --------------------------

% With var(x) = rms(x).^2 - mean(x).^2 we can compute the RMS value in dB with 
% the following formula:
% y = 20*log10( sqrt(var(x))*10^5 )
y = 10*log10( var(x)*10^10 );
