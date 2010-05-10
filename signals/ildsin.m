function outsig = ildsin(f,ild,fs)
%ILDSIN Generate a sinusoid with a interaural level difference (ILD)
%   Usage: outsig = itdsin(f,ild,fs)
%
%   Input parameters:
%       f       - carrier frequency of the sinusoid (Hz)
%       ild     - interaural level difference of the right signal, this can be
%                 positive or negative (dB)
%       fs      - sampling rate (Hz)
%
%   Output parameters:
%       outsig  - two channel 1 s long sinusoid (scaled to max(outsig)==1-eps)
%
%   ILDSIN(f,ild,fs) generates a sinusoid with a interaural level difference
%   of ild and a frequency of f.
%
%R moore2003introduction
%

% AUTHOR: Hagen Wierstorf

% ------ Checking of input parameters ---------

error(nargchk(3,3,nargin));

if ~isnumeric(f) || ~isscalar(f) || f<0
    error('%s: f must be a positive scalar.',upper(mfilename));
end

if ~isnumeric(ild) || ~isscalar(ild)
    error('%s: itd must be a scalar.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs must be a positive scalar!',upper(mfilename));
end


% ------ Computation --------------------------

% Create a one second time 
t = (1:fs)/fs;
% Left signal
sigl = sin(2*pi*f.*t);
% Right signal with level difference of ILD
sigr = gaindb(sigl,ild);
% Combine left and right channel
outsig = [sigl' sigr'];
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);
