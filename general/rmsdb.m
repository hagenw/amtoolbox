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
%   The level in dB is calculated by
%
%C 	level(x) = 20*log10( sqrt( rms(x).^2 - mean(x).^2 ) )
%
%   See also: rmsdb

%   AUTHOR : Hagen Wierstorf
  
% ------ Checking of input parameters ---------

error(nargchk(1,1,nargin));

if ~isnumeric(x) || ~isvector(x)
  error('RMS: Input must be a vector.');
end;

% ------ Computation --------------------------

% With var(x) = rms(x).^2 - mean(x).^2 we can compute the RMS value in dB with 
% the following formula:
y = 10*log10( var(x) );
