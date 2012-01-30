function [bmld_out] = culling2007bmld(coherence,phase_target,phase_int,fc)
%CULLING2007BLMD  Binaural Masking Level Difference 
%
%   Input parameters:
%      coherence    : (0->1) maximum of the interaural cross-correlation
%      phase_target : (0->2*pi) interaural phase difference of target signal 
%      phase_int    : (0->2*pi) interaural phase difference of the interfering source
%      fc           : Centre frequency of the frequency channel
%
%   Output parameters:
%      bmld_output  : improvement in predicted signal threshold when 
%                     *phase_target* and *phase_int* differ compared to 
%                     when they are the same.
%
%   `culling2007bmld(coherence,phase_target,phase_int,fc)` calculates the
%   binaural masking level differencefor a signal in broadband noise. 
%   The input noise coherence and phase must be pre-calculated for the 
%   frequency channel bearng the signal. See |jelfs2001|_ for an example on
%   how to calculate these.
% 
%   See also: jelfs2011
% 
%   References:  culling2007evidence
  
k = (1 + 0.25^2) * exp((2*pi*fc)^2 * 0.000105^2);
bmld_out = 10 * log10 ((k - cos(phase_target-phase_int))/(k - coherence));

bmld_out = max(mbld_out,0);

