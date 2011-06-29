function [outsig, fc] = breebaart2001preproc(insig, fs, varargin);
%BREEBAART2001PREPROC   Auditory model from Breebaart et. al. 2001
%   Usage: [outsig, fc] = breebaart2001preproc(insig,fs);
%          [outsig, fc] = breebaart2001preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
%   BREEBAART2001PREPROC(insig,fs) computes the internal representation of
%   the signal insig sampled with a frequency of fs Hz as described in
%   Breebaart (2001).
%  
%   [outsig,fc]=BREEBAART2001PREPROC(...) additionally returns the center
%   frequencies of the filter bank.
%
%   The following parameters may be passed at the end of the line of
%   input arguments:
%
%-     'flow',flow - Set the lowest frequency in the filterbank to
%                    flow. Default value is 80 Hz.
%
%-     'fhigh',fhigh - Set the highest frequency in the filterbank to
%                    fhigh. Default value is 8000 Hz.
%
%-     'basef',basef - Ensure that the frequency basef is a center frequency
%                    in the filterbank. The default value of [] means
%                    no default.
%
%   The model assumes than a pure tone input signal with an RMS value of 1
%   corresponds to an acoustic signal of 100 db SPL.
%  
%   The breebaart2001 model consists of the following stages:
%   
%     * a gammatone filter bank with 1-erb spaced filters.
%
%     * an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 770 Hz.
%
%     * an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%R  breebaart2001binaural

%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ------------

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.keyvals.flow=80;
definput.keyvals.fhigh=8000;
definput.keyvals.basef=[];

[flags,keyvals,flow,fhigh,basef]  = ltfatarghelper({'flow', 'fhigh', ...
                    'basef'},definput,varargin);

% ------ do the computation -------------------------

% find the center frequencies used in the filterbank, 1 ERB spacing
fc = erbspacebw(flow, fhigh, 1, basef);

% Calculate filter coefficients for the gammatone filter bank.
[gt_b, gt_a]=gammatone(fc, fs, 'complex');

% Apply the Gammatone filterbank
outsig = 2*real(filterbankz(gt_b,gt_a,insig));

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'breebaart');

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'breebaart');



