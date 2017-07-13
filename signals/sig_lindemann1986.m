function outsig = sig_lindemann1986(fc,mf,fs)
%sig_lindemann1986 Generate a binaural modulated sinus for the Lindemann (1986) model
%   Usage: outsig = sig_lindemann1986(fc,mf,fs)
%
%   Input parameters:
%       fc  : carrier frequency of the sinus (Hz)
%       mf  : binaural modulation frequency (Hz)
%       fs  : sampling rate (Hz)
%
%   Output parameters:
%       outsig  : $fs \times 2$ sinusoid signal
%
%   `sig_lindemann1986(fc,mf,fs)` generates an binaural modulated sinusoid with a
%   carrier frequency of *f* and a frequency moving around the two ears of
%   *mf*.
%
%   See also: lindemann1986 demo_lindemann1986

% AUTHOR: 13.04.2010 Hagen Wierstorf
%         27.01.2012 Peter Soendergaard
%         09.07.2017: Piotr Majdak

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

