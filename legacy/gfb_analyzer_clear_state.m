function analyzer = gfb_analyzer_clear_state(analyzer)
%GFB_ANALYZER_CLEAR_STATE  Reset filter states
%   Usage: analyzer = gfb_analyzer_clear_state(analyzer)
%
%   `analyzer=gfb_analyzer_clear_state(analyzer)` resets the filter states
%   to zeros
  
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006, Feb 2007

warning('Warning: GFB_ANALYZER_CLEAR_STATE will be removed in a future release. Use HOHMANN2002CLEARSTATE instead. ');

for band = 1:length(analyzer.center_frequencies_hz)
  analyzer.filters(1, band) = ...
      gfb_filter_clear_state(analyzer.filters(1, band));
end
