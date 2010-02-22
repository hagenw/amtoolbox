function outsig = bmsin(f,mf,fs)
%AMSIN Generate a binaural modulated sinus
%   Usage: outsig = bmsin(f,mf,fs)
%
%   Input parameters:
%       f       - carrier frequency of the sinus
%       mf      - binaural modulation frequency
%       fs      - sampling rate
%
%   Output parameters:
%       outsig  - two channel 1 s long sinusoid
%
%   BMSIN(f,mf,fs) generates an binaural modulated sinusoid with a
%   freuqnecy of f and a frequency moving around the two ears of mf.

% Create a one second time 
t = 1/fs:1/fs:1;
% Left signal
sigl = sin(2*pi*f.*t);
% Right signal with amplitude modulation
sigr = sin(2*pi*f.*t + sin(2*pi*mf.*t));
outsig = [sigl' sigr'];
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);