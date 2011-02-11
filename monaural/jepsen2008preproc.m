function [outsig, fc, mfc] = jepsen2008preproc(insig, fs, varargin);
%JEPSEN2008PREPROC   Auditory model from Jepsen et. al. 2008
%   Usage: [outsig, fc] = jepsen2008preproc(insig,fs);
%          [outsig, fc] = jepsen2008preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
%   JEPSEN2008PREPROC(insig,fs) computes the internal representation of the signal insig
%   sampled with a frequency of fs Hz as described in Jepsen, Ewert and
%   Dau (2008).
%  
%   [outsig,fc]=JEPSEN2008(...) additionally returns the center frequencies of
%   the filter bank.
%
%   The model assumes than a pure tone input signal with an RMS value of 1
%   corresponds to an acoustic signal of 100 db SPL.
%  
%   The Jepsen2008 model consists of the following stages:
% 
%     * A heaphone filter to simulate the effect of a standard set of
%     headphones
%
%     * A middle ear filter to simulate the effect of the middle ear, and
%     to convert to stapes movement.
%
%     * DRNL - Dual resonance non-linear filterbank
%
%     * an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz.
%
%     * An expansion stage
%
%     * an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%     * a modulation filterbank
%
%R  jepsen2008mapI jepsen2008mapII

%   AUTHOR : Torsten Dau, Morten Løve Jepsen, Peter L. Soendergaard
  
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

definput.import={'auditoryfilterbank'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig, fs, 'argimport',flags,keyvals);

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'dau');

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'dau');

% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);




