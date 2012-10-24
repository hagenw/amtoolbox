function [cqmag,fc,cqmaghr,fvec] = cqdft( insig,varargin )
%CQDFT FFT-based filter bank with constant relative bandwidth according
%   Usage:  [cqmag] = cqdft( insig )
%           [cqmag,fc,cqmaghr,fvec] = cqdft( insig,fs,flow,fhigh,bw )
%
%   Input parameters:
%      insig   : Impulse response or complex spectrum
%      fs      : Sampling rate, default is 48kHz.
%      flow    : Lowest frequency, minimum: 0.5kHz, default is 2kHz
%      fhigh   : Highest frequency, default is, default is 16kHz  
%      bw      : bandwidth, possible values 3,6,9,12, default is 6.
%
%   Output parameters:
%      cqmag   : mean magnitudes of CQ-bands in dB
%      fc      : center frequencies of bands (geo. mean of corners)
%      cqmaghr : same as cqmag but for all freq. bins (high resolution)
%      fvec    : freq. vector according to FFT-resolution
%
%   `cqdft(insig)` approximates a constant-Q filter bank by averaging the
%   magnitude bins of a DFT. `cqdft` results in 'bw' dB-magnitudes per octave.
%
%   References: langendijk2002contribution
    
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute

definput.keyvals.fs=48000;
definput.keyvals.flow=2000;
definput.keyvals.fhigh=16000;
definput.keyvals.bw=6;

[flags,kv]  = ltfatarghelper({'fs','flow','fhigh','bw'},definput,varargin);


% input signal given in time or frequency domain?
if isreal(insig)    % -> TD
    nfft = 2^12;%max(2^12,size(insig,1));
    y = abs(fft(insig,nfft));
else                % -> FD
    y = abs(insig);
    nfft = size(insig,1);
end
fvec = 0:kv.fs/nfft:kv.fs-kv.fs/nfft;

octs = log2(kv.fhigh/kv.flow);    % # of octaves
jj = 0:octs*kv.bw; 
n = round(2.^((jj)/kv.bw)*kv.flow/kv.fs*nfft);       % startbins
fc = zeros(length(jj)-1,1);                       % center frequencies
cqmag = zeros(length(jj)-1,size(y,2),size(y,3));  % mean magnitudes of CQ-bands
cqmaghr = zeros(size(y));                         % same but for all freq. bins (high resolution)
for ind = jj(1)+1:jj(end)
    nj = n(ind+1)-n(ind);
    idn = n(ind):n(ind+1)-1;
    fc(ind) = geomean(fvec(n(ind:ind+1)));
    cqmag(ind,:,:) = sqrt(1/(nj)*sum(y(idn,:,:).^2,1));
    cqmaghr(idn,:,:) = repmat(cqmag(ind,:,:),[length(idn),1,1]);
end

cqmag = db(cqmag);
cqmaghr(nfft/2+2:end,:,:) = cqmaghr(nfft/2:-1:2,:,:);
cqmaghr = db(cqmaghr);

end

