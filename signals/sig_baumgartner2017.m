function Obj = sig_baumgartner2017(Obj,M,flow,fhigh)
% sig_baumgartner2017 - flattens magnitude spectra of HRTFs.
%
%   Usage:  Obj = sig_baumgartner2017(Obj,M,flow,fhigh)
%
%   Input parameters:
%     Obj   : reference SOFA object
%     M     : spectral salience. 1 refers to reference, 0 to flat, -1 to
%             flipped spectral magnitude
%     flow  : lower cut-off frequency
%     fhigh : higher cut-off frequency
%
%   Output parameters:
%     Obj   : modified SOFA object

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

Nfft = 2^10;

ref = shiftdim(Obj.Data.IR,2);

hrtf = fftreal(ref,Nfft);
freq = 0:Obj.Data.SamplingRate/Nfft:Obj.Data.SamplingRate/2; % frequency vector

% following proecessing only done for dedicated frequency range
idf = freq >= flow & freq <= fhigh; % indices for dedicated frequency range
Nf = sum(idf); % # frequency bins
mag = db(abs(hrtf(idf,:,:))); % HRTF magnitudes in dB
idwf = idf(:) | circshift(idf(:),[1,0]); % include one neighbouring position for evaluation of frequency weighting
wf = diff(freqtoerb(freq(idwf))); % frequency weighting according to differentiated ERB scale
wf = repmat( wf(:)/sum(wf) ,1,size(mag,2),size(mag,3));
meanmag = repmat(sum(wf.*mag,1),Nf,1,1);
varmag = mag - meanmag;

ph = angle(hrtf);
modmag = db(abs(hrtf));
modmag(idf,:,:) = meanmag + M*varmag;
mod = ifftreal(10.^(modmag/20).*exp(1i*ph),Nfft);

Obj.Data.IR = shiftdim(mod,1);

end