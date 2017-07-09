function outsig = sig_ildsin(fc,ild,fs)
%sig_ildsin Sinusoid with a interaural level difference (ILD)
%   Usage: outsig = sig_itdsin(fc,ild,fs)
%
%   Input parameters:
%       fc      : carrier frequency of the sinusoid (Hz)
%       ild     : ILD of the right signal, positive or negative (dB)
%       fs      : sampling rate (Hz)
%
%   Output parameters:
%       outsig  : two channel 1 s long sinusoid 
%
%   `sig_ildsin(fc,ild,fs)` generates a sinusoid with a interaural level difference
%   of *ild* and a frequency of *fc*.
%
%   The output is scaled to have a maximum value of 1-eps.
%
%   References: moore2003introduction

% AUTHOR: Hagen Wierstorf

% ------ Checking of input parameters ---------

error(nargchk(3,3,nargin));

if ~isnumeric(fc) || ~isscalar(fc) || fc<0
    error('%s: fc must be a positive scalar.',upper(mfilename));
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
sigl = sin(2*pi*fc.*t);
% Right signal with level difference of ILD
sigr = gaindb(sigl,ild);
% Combine left and right channel
outsig = [sigl' sigr'];
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);

