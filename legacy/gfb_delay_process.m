function [outsig, delay] = gfb_delay_process(delay, insig)
%GFB_DELAY_PROCESS  Filterbank delay processing
%   Usage: [outsig, delay] = gfb_delay_process(delay, insig)
%
%   Input parameters:
%     delay  : A `gfb_delay` structure created from |gfb_delay_new|. The delay
%              will be returned with updated delayline states as the second
%              return parameter
%     insig  : A complex matrix containing the signal to delay.  Each row
%              corresponds to a filterbank band
%
%   Output parameters:
%     outsig : A real matrix containing the delay's output
%
%   `gfb_delay_process(delay, insig)` will delay each band (row) of the
%   input data *insig* by a band-dependent amount of samples, then
%   multiplied with a band-dependent complex constant.  Finally, the real
%   part of this product will be returned.
%
%   See also: gfb_delay_new

% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

warning('Warning: GFB_DELAY_PROCESS will be removed in a future release. Use hohmann2002_process instead. ');

[number_of_bands, number_of_samples] = size(insig);
if (number_of_bands ~= length(delay.delays_samples))
  error('input rows must match the number of bands');
end
outsig = zeros(number_of_bands, number_of_samples);
for band = [1:number_of_bands]
  if (delay.delays_samples(band) == 0)
    outsig(band,:) = ...
        real(insig(band,:) * delay.phase_factors(band));
  else
    tmp_out = [delay.memory(band,1:delay.delays_samples(band)), ...
               real(insig(band,:) * delay.phase_factors(band))];
    delay.memory(band,1:delay.delays_samples(band)) = ...
        tmp_out(number_of_samples+1:length(tmp_out));
    outsig(band,:) = tmp_out(1:number_of_samples);
  end
end

