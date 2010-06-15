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
%   DAU96PREPROC(insig,fs) computes the internal representation of the signal insig
%   sampled with a frequency of fs Hz as described in Dau, Puschel and
%   Kohlrausch (1996a).
%  
%   [outsig,fc]=DAU96(...) additionally returns the center frequencies of
%   the filter bank.
%
%   The following parameters may be passed at the end of the line of
%   input arguments:
%
%      'flow',flow - Set the lowest frequency in the filterbank to
%                    flow. Default value is 80 Hz.
%
%      'fhigh',fhigh - Set the highest frequency in the filterbank to
%                    fhigh. Default value is 8000 Hz.
%
%      'basef',basef - Ensure that the frequency basef is a center frequency
%                    in the filterbank. The default value of [] means
%                    no default.
% 
%      'subfs',subfs - Apply a final downsampling of the subband signals
%                    to subfs Hz to avoid excessive data. The default value
%                    of [] means no downsampling.
%
%   The model assumes than a pure tone input signal with an RMS value of 1
%   corresponds to an acoustic signal of 100 db SPL.
%  
%   The dau96 model consists of the following stages:
%   
%     * a gammatone filter bank with 1-erb spaced filtes.
%
%     * an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz.
%
%     * an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%     * a modulation low pass filter liming modulations to below 50 Hz.
%
%   The model implemented in this file is not identical to the model
%   published in Dau et. al. (1996a). An overshoot limit has been added to
%   the adaptation stage to fix a problem where abrupt changes in the
%   input signal would cause unnaturally big responses. This is described
%   in Dau et. al. (1997a).
%
%R  dau1996qmeI dau1996qmeII dau1997mapI

%   AUTHOR : Torsten Dau, Morten, Løve Jepsen, Peter L. Soendergaard
  
% ------ Checking of input parameters ------------

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(inoutsig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.keyvals.flow=80;
definput.keyvals.fhigh=8000;
definput.keyvals.basef=-1;
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({},defnopos,varargin);

% ------ do the computation -------------------------

% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(flags.flow, flags.fhigh, 1, flags.basef);

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
if ~isempty(flags.subfs)
  inoutsig = fftresample(inoutsig,round(length(inoutsig)/fs*subfs));
end;


