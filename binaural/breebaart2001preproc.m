function [ei_map, fc] = breebaart2001preproc(insig, fs, tau, ild, varargin);
%BREEBAART2001PREPROC   Auditory model from Breebaart et. al. 2001
%   Usage: [outsig, fc] = breebaart2001preproc(insig,fs);
%          [outsig, fc] = breebaart2001preproc(insig,fs,...);
%
%   Input parameters:
%        insig  : input acoustic signal.
%        fs     : sampling rate.
%        tau    : characteristic delay in seconds (positive: left is leading)
%        ild    : characteristic ILD in dB (positive: left is louder)
%  
%   BREEBAART2001PREPROC(insig,fs,tau,ild) computes the EI-cell representation of
%   the signal insig sampled with a frequency of fs Hz as described in
%   Breebaart (2001). The parameters tau and ild define the sensitivity of the EI-cell.
%
%   The input must have dimensions time x left/right channel x signal no.
%
%   The output has dimensions time x frequency x signal no. 
%  
%   [outsig,fc]=BREEBAART2001PREPROC(...) additionally returns the center
%   frequencies of the filter bank.
%  
%   The Breebaart 2001 model consists of the following stages:
%   
%     1) a gammatone filter bank with 1-erb spaced filters.
%
%     2) an envelope extraction stage done by half-wave rectification
%        followed by low-pass filtering to 770 Hz.
%
%     3) an adaptation stage modelling nerve adaptation by a cascade of 5
%        loops.
%
%     4) An excitation-inhibition (EI) cell model.
%
%   Parameters for AUDITORYFILTERBANK, IHCENVELOPE, ADAPTLOOP and EICELL can be
%   passed at the end of the line of input arguments.
%
%   See also: eicell, auditoryfilterbank, ihcenvelope, adaptloop

%R  breebaart2001binaural

%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ------------

if nargin<4
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import = {'auditoryfilterbank','ihcenvelope','adaptloop','eicell'};
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
outsig = 2*real(ufilterbankz(gt_b,gt_a,insig,1));

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'ihc_breebaart');

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'adt_breebaart');

[siglen,nfreqchannels,naudiochannels,nsignals] = size(outsig);

ei_map = zeros(siglen, nfreqchannels, nsignals);
for k=1:nsignals
  for g=1:nfreqchannels
    ei_map(:,g,k) = eicell(squeeze(outsig(:,g,:,k)),fs,tau,ild);
  end
end


