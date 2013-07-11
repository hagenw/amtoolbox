function [outsig, fc, mfc] = dau1997preproc(insig, fs, varargin);
%DAU1997PREPROC   Auditory model from Dau et. al. 1997
%   Usage: [outsig, fc] = dau1997preproc(insig,fs);
%          [outsig, fc] = dau1997preproc(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
%   `dau1997preproc(insig,fs)` computes the internal representation of the
%   signal *insig* sampled with a frequency of *fs* Hz as described in Dau,
%   Puschel and Kohlrausch (1997a).
%  
%   `[outsig,fc,mfc]=dau1997preproc(...)` additionally returns the center
%   frequencies of the filter bank and the center frequencies of the
%   modulation filterbank.
%  
%   The Dau 1997 model consists of the following stages:
%   
%   1) a gammatone filter bank with 1-erb spaced filtes.
%
%   2) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 1000 Hz.
%
%   3) an adaptation stage modelling nerve adaptation by a cascade of 5
%      loops.
%
%   4) a modulation filterbank
%
%   Any of the optinal parameters for |auditoryfilterbank|,
%   |ihcenvelope| and |adaptloop| may be optionally specified for this
%   function. They will be passed to the corresponding functions.
%
%   See also: auditoryfilterbank, ihcenvelope, adaptloop, modfilterbank

%   References: dau1997mapI dau1997mapII

%   AUTHOR : Torsten Dau, Morten Løve Jepsen, Peter L. Søndergaard
  
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

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop'};
definput.importdefaults={'ihc_dau','adt_dau'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(insig,fs,'argimport',flags,keyvals);

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,keyvals);

% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);

