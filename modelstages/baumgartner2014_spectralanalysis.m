function varargout = baumgartner2014_spectralanalysis(sig,varargin)
%baumgartner2014_spectralanalysis - Approximation of spectral analysis by auditory periphery
%   Usage:     mp = baumgartner2014_spectralanalysis(sig)
%
%   Input parameters:
%     sig     : incoming time-domain signal
%
%   Output parameters:
%     mp      : spectral magintude profile
%     fc      : center frequencies of auditory filters
%
%   `baumgartner2014_spectralanalysis(...)` computes temporally integrated
%   spectral magnitude profiles.
%
%   `baumgartner2014_spectralanalysis` accepts the following optional parameters:
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    *flow*. Default value is 700 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    *fhigh*. Default value is 18000 Hz.
%
%     'fs',fs        Define the sampling rate of the impulse responses. 
%                    Default value is 48000 Hz.
%
%     'space',sp     Set spacing of auditory filter bands (i.e., distance 
%                    between neighbouring bands) to *sp* in number of
%                    equivalent rectangular bandwidths (ERBs). 
%                    Default value is 1 ERB.
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

definput.import={'baumgartner2014'};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% Spectral Analysis, Eq.(2)

if kv.space == 1 % Standard spacing of 1 ERB
  [mp,fc] = auditoryfilterbank(sig(:,:),kv.fs,'flow',kv.flow,'fhigh',kv.fhigh);
else
  fc = audspacebw(kv.flow,kv.fhigh,kv.space,'erb');
  [bgt,agt] = gammatone(fc,kv.fs,'complex');
  mp = 2*real(ufilterbankz(bgt,agt,sig(:,:)));  % channel (3rd) dimension resolved
end
Nfc = length(fc);   % # of bands

% Set back the channel dimension
mp = reshape(mp,[size(sig,1),Nfc,size(sig,2),size(sig,3)]);

% Averaging over time (RMS)
mp = 20*log10(squeeze(rms(mp)));      % in dB

if size(mp,2) ~= size(sig,2) % retreive polar dimension if squeezed out
    mp = reshape(mp,[size(mp,1),size(sig,2),size(sig,3)]);
end

varargout{1} = mp;
if nargout == 2
  varargout{2} = fc;
end

end