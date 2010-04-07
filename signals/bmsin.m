function outsig = bmsin(f,mf,fs)
%BMSIN Generate a binaural modulated sinus
%   Usage: outsig = bmsin(f,mf,fs)
%
%   Input parameters:
%       f       - carrier frequency of the sinus (Hz)
%       mf      - binaural modulation frequency (Hz)
%       fs      - sampling rate (Hz)
%
%   Output parameters:
%       outsig  - two channel 1 s long sinusoid
%
%   BMSIN(f,mf,fs) generates an binaural modulated sinusoid with a
%   carrier frequency of f and a frequency moving around the two ears of mf.
%

% AUTHOR: Hagen Wierstorf

% ------ Checking of input parameters ---------

error(nargchk(3,3,nargin));

if ~isnumeric(f) || ~isscalar(f) || f<0
    error('%s: f must be a positive scalar.',upper(mfilename));
end

if ~isnumeric(mf) || ~isscalar(mf) || mf<=0
    error('%s: mf must be a positive scalar.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs must be a positive scalar!',upper(mfilename));
end


% ------ Computation --------------------------
% Create a one second time 
t = (1:fs)/fs;
% Left signal
sigl = sin(2*pi*f.*t);
% Right signal with amplitude modulation
sigr = sin(2*pi*f.*t + sin(2*pi*mf.*t));
outsig = [sigl' sigr'];
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);
