function center_frequencies_hz = gfb_center_frequencies(filters_per_ERBaud,flow,basef,fhigh)
%GFB_CENTER_FREQUENCIES  Obsolete, use erbspacebw instead
%   Usage: frequencies_hz = gfb_center_frequencies(frequencies_per_ERBaud,flow,basef,fhigh);
% 
%   The call::
%
%     gfb_center_frequencies(filters_per_ERBaud,flow,basef,fhigh)
%
%   can be replaced by::
%
%     erbspacebw(flow,fhigh,1/filters_per_ERBaud,basef);
%
  
if (nargin < 4)
  fhigh = basef;
end

% Calculate the values of the parameter frequencies on the ERBscale:
lower_cutoff_frequency_erb     = freqtoerb(flow);
specified_center_frequency_erb = freqtoerb(basef);
upper_cutoff_frequency_erb     = freqtoerb(fhigh);


% The center frequencies of the individual filters are equally
% distributed on the ERBscale.  Distance between adjacent filters'
% center frequencies is 1/filters_per_ERBaud.
% First, we compute how many filters are to be placed at center
% frequencies below the base frequency:
erbs_below_base_frequency = ...
    specified_center_frequency_erb - lower_cutoff_frequency_erb;
num_of_filters_below_base_freq = ...
    floor(erbs_below_base_frequency * filters_per_ERBaud);

% Knowing this number of filters with center frequencies below the
% base frequency, we can easily compute the center frequency of the
% gammatone filter with the lowest center frequency:
start_frequency_erb = ...
    specified_center_frequency_erb - ...
    num_of_filters_below_base_freq / filters_per_ERBaud;

% Now we create a vector of the equally distributed ERBscale center
% frequency values:
center_frequencies_erb = ...
    [start_frequency_erb:(1/filters_per_ERBaud):upper_cutoff_frequency_erb];
center_frequencies_hz = erbtofreq(center_frequencies_erb);
