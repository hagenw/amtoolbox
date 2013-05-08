function [bmld_out] = culling2005bmld(coherence,phase_target,phase_int,fc)
%CULLING2005BMLD  Binaural Masking Level Difference 
%
%   Input parameters:
%      coherence    : Maximum of the interaural cross-correlation, the
%                     value should be between 0 and 1.
%      phase_target : Interaural phase difference of target signal,
%                     between 0 and $2\pi$
%      phase_int    : Interaural phase difference of the interfering
%                     source, between 0 and $2\pi$
%      fc           : Centre frequency of the frequency channel
%
%   Output parameters:
%      bmld_output  : improvement in predicted signal threshold when 
%                     *phase_target* and *phase_int* differ compared to 
%                     when they are the same.
%
%   `culling2005bmld(coherence,phase_target,phase_int,fc)` calculates the
%   binaural masking level differencefor a signal in broadband noise. 
%   The input noise coherence and phase must be pre-calculated for the 
%   frequency channel bearng the signal. See |jelfs2011| for an example on
%   how to calculate these.
% 
%   See also: jelfs2011
% 
%   References:  culling2005erratum culling2004role
  
k = (1 + 0.25^2) * exp((2*pi*fc)^2 * 0.000105^2);
bmld_out = 10 * log10 ((k - cos(phase_target-phase_int))/(k - coherence));

bmld_out = max(bmld_out,0);

