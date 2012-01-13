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
%     1) A heaphone filter to simulate the effect of a standard set of
%     headphones
%
%     2) A middle ear filter to simulate the effect of the middle ear, and
%     to convert to stapes movement.
%
%     3) DRNL - Dual resonance non-linear filterbank
%
%     4) an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 1000 Hz.
%
%     5) An expansion stage
%
%     6) an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%     7) a modulation filterbank
%
%   References:jepsen2008cmh

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

definput.import={'drnl'};
definput.importdefaults={'jepsen2008'};
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

%% Headphone filter
hp_fir = headphonefilter(fs);
outsig = filter(hp_fir,1,insig);

%% DRNL and compensation for middle-ear
[outsig, fc] = drnl(outsig, fs, 'argimport',flags,keyvals);
outsig = gaindb(outsig,50);

%% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'ihc_dau');

%% Expansion stage
outsig = outsig.^2;

%% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'adt_dau');

%% Modulation filterbank
[outsig,mfc] = modfilterbank(outsig,fs,fc);

%OLDFORMAT
