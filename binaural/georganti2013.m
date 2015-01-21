function difSTD=georganti2013(signal,P)
%GEORGANTI2013 Binaural spectral-magnitude difference standard deviation according to Georganti et al., (2013)
%   Usage: difSTD = georganti2013(signal, P)
%
%   Input parameters:
%       signal : binaural input signal
%
%       P.fs: sampling rate in Hz
%
%       P.timeFr: Frame size in seconds
%
%       P.fmin: lower frequency (Hz) for the BSDM STD calculation
%
%       P.fmax: upper frequency (Hz) for the BSDM STD calculation
%
%   Output parameters:
%       difSTD: Binaural spectral-magnitude difference standard deviation (dB)
%
%   See also: exp_georganti2013
%
%   References: georganti2013 georganti2013application
%

%
%   Developed with Matlab 7.1.0.584 (R2010b) v.0.1
%   Updated with Matlab R2011b (7.13.0.564) v.1.0
%
%   Please send bug reports to:
%     Eleftheria Georganti
%     Postdoctoral Researcher
%     Experimental Audiology, ENT
%     University Hospital of Zurich/University of Zurich
%     Zurich, Switzerland
%     eleftheria.georganti@uzh.ch
%
%



if ~exist('P','var'), P=[]; end
  
if ~isfield(P,'fmin'), P.fmin=20; end % lower frequency in Hz - default value
if ~isfield(P,'fmax'), P.fmax=23000; end  % upper frequency in Hz - default value
if ~isfield(P,'fs'), P.fs=44100; end
if ~isfield(P,'timeFr'), P.timeFr=1; end

P.sampleFr = P.fs * P.timeFr; % Frame size in samples
if ~isfield(P,'hop'), P.hop = P.sampleFr/2; end % Overlap
P.nFFT = P.sampleFr; % FFT points

P.freq = (P.fs/P.nFFT)*(0:(P.nFFT/2-1)); % Frequency index


fmin_id = min(find((P.freq>P.fmin)));
fmax_id = min(find((P.freq>P.fmax)));

difSTD=zeros(1,length(1:P.hop:length(signal)-P.hop));
idx = 1;

for kk = 1:P.hop:length(signal)-P.hop

    % Calculate magnitude spectrums in dB of the left & right signals
    leftFFT  = 20*log10(abs(fft(signal(kk:kk+P.hop-1,1))));
    rightFFT = 20*log10(abs(fft(signal(kk:kk+P.hop-1,2))));

    % Subtract the magnitude spectrums
    specDIF  = leftFFT(1:end/2)-rightFFT(1:end/2);

    % Calculate the differential standard deviation for the
    % frequency range of interest
    difSTD(1,idx) = std(specDIF(fmin_id:fmax_id));

    idx = idx+1;

end
    

