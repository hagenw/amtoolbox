function [inoutsig, fc] = dau96preproc(inoutsig, fs, flow, fhigh,subfs,basef)
%DAU96PREPROC   Auditory model from Dau et. al. 1996.
%   Usage: [outsig, fc] = dau96preproc(insig,fs,flow,fhigh);
%          [outsig, fc] = dau96preproc(insig,fs,flow,fhigh,subfs);
%          [outsig, fc] = dau96preproc(insig,fs,flow,fhigh,subfs,basef);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%     flow   : lowest filter center frequency. 
%     fhigh  : highest filter center frequency.
%     basef  : Always include this frequency in the filter bank (optional).
%  
%   DAU96PREPROC(insig,fs,flow,fhigh) computes the internal representation of
%   the signal insig sampled with a frequency of fs Hz as described in Dau,
%   Puschel and Kohlrausch (1996a). The parameters flow and fhigh
%   determine the lowest and highest frequencies in the filter bank. The
%   filters are spaced 1 Erb apart. 
%
%   DAU96PREPROC(insig,fs,flow,fhigh,basef) does the same, but ensures that the
%   frequency basef is one of the frequencies in the filter bank.
%
%   DAU96PREPROC(insig,fs) chooses the lowest frequency to be 80 and the highest
%   to be 8000.
%
%   [outsig,fc]=DAU96PREPROC(...) additionally returns the center frequencies of
%   the filter bank.
%
%   The model assumes than an input signal with an RMS value of 1
%   corresponds to an acoustic signal of 100 db SPL.
%  
%   The dau96 model consists of the following stages:
%   
%     * a gammatone filter bank with 1-erb spaced filtes.
%
%     * an envelope extraction stage done by half-wave rectification
%          followed by low-pass filtering to 1000 Hz.
%
%     * an adaptation stage modelling nerve adaptation by a cascade of 5
%          loops.
%
%     * a modulation low pass filter liming modulations to below 50 Hz.
%-
%
%   The model implemented in this file is not identical to the model
%   published in Dau et. al. (1996a). An overshoot limit has been added to
%   the adaptation stage to fix a problem where abrupt changes in the
%   input signal would cause unnaturally big responses. This is described
%   in Dau et. al. (1997a).
%
%   See also: dau97, jepsen08
%
%R  dau1996qmeI dau1996qmeII dau1997mapI

%   AUTHOR : Torsten Dau, Morten, Løve Jepsen, Peter L. Soendergaard
  
% ------ Checking of input parameters ------------

error(nargchk(2,6,nargin));

if nargin<6
  basef=-1
end;

if nargin<5
  subfs=100;
end;

if nargin==2
  flow=80;
  fhigh=8000;
end;

if ~isnumeric(inoutsig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(flow) || ~isscalar(flow) || flow<0
  error('%s: flow must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) || fhigh<0
  error('%s: fhigh must be a non-negative scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;


% ------ do the computation -------------------------

% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(flow, fhigh, 1, basef);

% Calculate filter coefficients for the gammatone filter bank.
[gt_b, gt_a]=gammatone(fc, fs);

% Apply the Gammatone filterbank
inoutsig = 2*real(filterbank(gt_b,gt_a,inoutsig));

% 'haircell' envelope extraction
inoutsig = ihcenvelope(inoutsig,fs,'dau');

% non-linear adaptation loops
inoutsig = adaptloop(inoutsig,fs,'dau');

% Calculate filter coefficients for the 20 ms (approx.eq to 8 Hz) modulation
% lowpass filter.
mlp_a = exp(-(1/0.02)/fs);
mlp_b = 1 - mlp_a;
mlp_a = [1, -mlp_a];

% Apply the low-pass modulation filter.
inoutsig = filter(mlp_b,mlp_a,inoutsig);

% Apply final resampling to avoid excessive data
inoutsig = fftresample(inoutsig,round(length(inoutsig)/fs*subfs));

