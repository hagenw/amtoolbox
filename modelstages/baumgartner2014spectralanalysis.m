function varargout = baumgartner2014spectralanalysis(sig,kv,flags)
%BAUMGARTNER2014SPECTRALANALYSIS - Approximation of spectral analysis by auditory periphery
%   Usage:     mp = baumgartner2014spectralanalysis(sig)
%              [mp,fc] = baumgartner2014spectralanalysis(sig,kv,flags)
%
%   Input parameters:
%     sig     : incoming time-domain signal
%
%   Output parameters:
%     mp      : spectral magintude profile
%     fc      : center frequencies of auditory filters
%
%   `baumgartner2014spectralanalysis(...)` computes temporally integrated
%   spectral magnitude profiles.
%
%   `baumgartner2014spectralanalysis` accepts the following optional parameters:
%
%     'flags',flags  Transfer flags. If not defaults of baumgartner2014 are used.
%
%     'kv',kv        Transfer key-value pairs. If not defaults of baumgartner2014 are used.
%
%   References: baumgartner2014modeling

% AUTHOR: Robert Baumgartner

if not(exist('flags','var') || exist('kv','var'))
  definput.import={'baumgartner2014'};
  [flags,kv]=ltfatarghelper({},definput,{});
end

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