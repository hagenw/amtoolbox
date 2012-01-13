function outsig = bmsin(fc,mf,fs)
%BMSIN Generate a binaural modulated sinus
%   Usage: outsig = bmsin(fc,mf,fs)
%
%   Input parameters:
%       fc  - carrier frequency of the sinus (Hz)
%       mf  - binaural modulation frequency (Hz)
%       fs  - sampling rate (Hz)
%
%   Output parameters:
%       outsig  - fs x 2 sinusoid signal
%
%   BMSIN(fc,mf,fs) generates an binaural modulated sinusoid with a
%   carrier frequency of f and a frequency moving around the two ears of mf.

% AUTHOR: Hagen Wierstorf

% ------ Checking of input parameters ---------

error(nargchk(3,3,nargin));

if ~isnumeric(fc) || ~isscalar(fc) || fc<0
    error('%s: f has to be a positive scalar.',upper(mfilename));
end

if ~isnumeric(mf) || ~isscalar(mf) || mf<=0
    error('%s: mf has to be a positive scalar.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end


% ------ Computation --------------------------
% Create a one second time 
t = (1:fs)/fs;
% Left signal
sigl = sin(2*pi*fc.*t);
% Right signal with amplitude modulation
sigr = sin(2*pi*fc.*t + sin(2*pi*mf.*t));
outsig = [sigl' sigr'];
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);

%OLDFORMAT
