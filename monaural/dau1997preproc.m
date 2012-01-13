function [outsig, fc, mfc] = dau1997preproc(insig, fs, varargin);
%DAU1997PREPROC   Auditory model from Dau et. al. 1997.
%   Usage: [outsig, fc] = dau1997preproc(insig,fs);
%          [outsig, fc] = dau1997preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
%   DAU1997PREPROC(insig,fs) computes the internal representation of the signal insig
%   sampled with a frequency of fs Hz as described in Dau, Puschel and
%   Kohlrausch (1996a).
%  
%   [outsig,fc]=DAU1997(...) additionally returns the center frequencies of
%   the filter bank.
%
%   The model assumes than a pure tone input signal with an RMS value of 1
%   corresponds to an acoustic signal of 100 db SPL.
%  
%   The Dau1997 model consists of the following stages:
%   
%     1) a gammatone filter bank with 1-erb spaced filtes.
%
%     2) an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz.
%
%     3) an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%     4) a modulation filterbank
%
%   References:dau1997mapI dau1997mapII

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
outsig = ihcenvelope(outsig,fs,'ihc_dau');

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'adt_dau');

% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);





%OLDFORMAT
