function [bmld_out] = culling2007bmld(coherence,phase_target,phase_int,fc)
%CULLING2007BLMD  Binaural Masking Level Difference 
%
%   `culling2007bmld(coherence,phase_target,phase_int,fc)` calculates the
%   binaural masking level difference. XXX
% 
%   See also: culling2010
% 
%   References:  culling2007evidence
  
k = (1 + 0.25^2) * exp((2*pi*fc)^2 * 0.000105^2);
bmld_out = 10 * log10 ((k - cos(phase_target-phase_int))/(k - coherence));
if bmld_out < 0;
    bmld_out = 0;
end
return
