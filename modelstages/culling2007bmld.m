function [bmld_out] = culling2007bmld(coherence,phase_target,phase_int,fc)
%CULLING2007BLMD  Binaural Masking Level Difference 
%
%   Input parameters:
%      coherence    : XXX
%      phase_target : XXX
%      phase_int    : XXX
%      fc           : Centre frequency of XXX
%
%   Output parameters:
%      bmld_output  : XXX
%
%   `culling2007bmld(coherence,phase_target,phase_int,fc)` calculates the
%   binaural masking level difference. XXX
% 
%   See also: jelfs2011
% 
%   References:  culling2007evidence
  
k = (1 + 0.25^2) * exp((2*pi*fc)^2 * 0.000105^2);
bmld_out = 10 * log10 ((k - cos(phase_target-phase_int))/(k - coherence));

bmld_out = max(mbld_out,0);

