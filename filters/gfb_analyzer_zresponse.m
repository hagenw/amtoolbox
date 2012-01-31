function zresponse = gfb_analyzer_zresponse(analyzer, z)
%GFB_ANALYZER_ZRESPONSE  Frequency response
%   USAGE: zresponse = gfb_analyzer_zresponse(analyzer, z)
%
%   Input parameters:
%     analyzer : A gfb_analyzer struct as created by gfb_analyzer_new.
%     z        : A vector of z-plane frequencies where the frequency response
%                should be computed. `z = exp(2i*pi*f/fs)`
%
%   Output parameters:
%     zresponse : The complex frequency response of the filter
%
%   `gfb_analyzer_zresponse(analyzer,z)` computes the frequency response of
%   the gammatone filters in the filterbank at the frequencies *z*.
%   The frequency responses of the filters in the filterbank are stored
%   as columns in the output.

% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan & Nov 2006, Jan Feb 2007

number_of_bands = length(analyzer.center_frequencies_hz);
z = z(:);
zresponse = ones(length(z), number_of_bands);

for band = [1:number_of_bands]
  filter = analyzer.filters(band);
  zresponse(:,band) = gfb_filter_zresponse(filter, z);
end
