function [weightings] = siiweightings(fc)
%SIIWEIGHTINGS  Compute the SII weightings
%   Usage: [weightings] = siiweightings(fc)
%
%   ` siiweightings(fc)` computes the SII-weighting for the centre
%   frequencies given in *fc*.
  
weightings = zeros(length(fc),1);
bands = [0, 100, 200, 300, 400, 4400, 5300, 6400, 7700, 9500].';
weights = [0, 0.0103, 0.0261, 0.0419, 0.0577, 0.0460, 0.0343, 0.0226, 0.0110, 0].';
for n = 1:length(fc)
  if fc(n) >= 9500
    weightings(n) = 0;
  else
    ii = find(bands > fc(n));
    weightings(n) = weights((ii(1) - 1));
  end
end
weightings = weightings ./ sum(weightings);
