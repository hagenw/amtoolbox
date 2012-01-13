function [benefit, weighted_SNR, weighted_bmld] = culling2010(target,interferer,fs,varargin)
%CULLING2010  Binaural advantage for speech in reverberant conditions
%   Usage:  [benefit weighted_SNR weighted_bmld] = culling2010(target,interferer,fs)
%  
%   Input parameters:
%     target       : Binaural target stimuli
%     interfererer : Binaural interferer stimuli
%    
%   `culling2010(target,interferer,fs)` computes the increase in speech
%   intelligibility of the target when listening binaurally to the target
%   and interferer. The signals are assumed to be sampled at a sampling
%   frequency of *fs* Hz.
%
%   See also: culling2007bmld
% 
%   References:  culling2010mapping culling2007evidence
  
definput.flags.ears={'both','left','right'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Make sure that there is slightly more than 1 erb per channel, and get
% the gammatone filters.
nchannels=ceil(freqtoerb(fs/2));
fc=erbspace(0,fs/2,nchannels);
[b,a] = gammatone(fc,fs,'complex');

effective_SNR = zeros(nchannels,1);
bmld_prediction = zeros(nchannels,1);

targ_f = 2*real(filterbankz(b,a,target));
int_f  = 2*real(filterbankz(b,a,interferer));

for n = 1:nchannels
  % Calculate the effect of BMLD
  if flags.do_both
    % cross-correlate left and right signal in channel n for both the
    % target and the inteferer
    [phase_t, coher_t] = do_xcorr(targ_f(:,n,1),targ_f(:,n,2),fs,fc(n)); 
    [phase_i, coher_i] = do_xcorr( int_f(:,n,1), int_f(:,n,2),fs,fc(n)); 

    bmld_prediction(n) = culling2007bmld(coher_t,phase_t,phase_i,fc(n));
  end

  % Calculate the effect of better-ear SNR
  left_SNR  = sum(targ_f(:,n,1).^2) / sum(int_f(:,n,1).^2);
  right_SNR = sum(targ_f(:,n,2).^2) / sum(int_f(:,n,2).^2);

  if flags.do_both
    SNR = max(left_SNR,right_SNR);
  end;

  if flags.do_left
    SNR = left_SNR;
  end;
  
  if flags.do_right
    SNR = right_SNR;
  end
  
  % combination
  effective_SNR(n) = 10*log10(SNR);
end

% Calculate the SII weighting
weightings = siiweightings(fc);

if flags.do_both
  weighted_bmld = sum(bmld_prediction.*weightings);
else
  weighted_bmld = 0;
end

weighted_SNR = sum(effective_SNR.*weightings);

benefit = weighted_SNR + weighted_bmld;


% Helper function to do the cross-correlation, and extract the delay of
% the peak (output parameter 'phase' and the coherence at the peak.
function [phase coherence] = do_xcorr(left, right, fs, fc)
  % Compute maximum lag to be considered. This could be very large for a
  % centre frequency close to 0, so set a maximum.
  maxlag = min(length(left),round(fs/(fc*2)));
  
  % Compute the correlation
  iacc = xcorr(squeeze(left),squeeze(right),maxlag,'coeff');
  
  % Find the position of the large correlation coefficient.
  [coherence delay_samp] = max(iacc);
  phase = fc*2*pi*delay_samp/fs;