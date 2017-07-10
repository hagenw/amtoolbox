function out = sig_boyd2012(in,source,varargin)
%sig_boyd2012 - Stimulus from Boyd et al. (2012) 
%   Usage: out = sig_boyd2012(in,source,mix,azi)
%
%   Time-aligned mixture between individualized binaural stimulus and 
%   head-absent simulation (ITD only).
%
%   Input parameters:
%     in      : binaural input signal.
%     source	: monaural source signal.
%     mix     : mixing ratio ranging from 0 to 1. Default is 1 and means 
%               out is same as in.
%     lp      : low-pass cut off frequency. Default is NAN and means
%               broadband.
%     fs      : sampling rate in Hz (only required for low-pass filtering).
%               Default is 48 kHz.
%     azi     : azimuth (positive to the left).
%
%   Output parameters:
%     out     : binaural output signal mixture.
%
%   References: boyd2012

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

definput.keyvals.mix = 2;
definput.keyvals.lp = nan;
definput.keyvals.fs = 48e3;
definput.keyvals.azi = -30;

[flags,kv]  = ltfatarghelper({'mix','lp','fs','azi'},definput,varargin);

% create time-aligned head-absent simulation
[scor,lag] = xcorr(in(:,1),source(:));
[~,iL] = max(abs(scor));
iL = lag(iL);
[~,iR] = max(abs(xcorr(in(:,2),source(:))));
iR = lag(iR);
ha = zeros(length(in),2);
ha(iL+(1:length(source)),1) = source;
ha(iR+(1:length(source)),2) = source;
if kv.azi > 0
  ha = fliplr(ha);
end
SPL = mean(dbspl(in));
ha = setdbspl(ha,SPL);
   
% mixing
out = kv.mix*in + (1-kv.mix)*ha;

% low-pass filtering
if not(isnan(kv.lp) || isempty(kv.lp))
  [b,a]=butter(10,2*kv.lp/kv.fs,'low');
  out = filter(b,a,out);
end

end