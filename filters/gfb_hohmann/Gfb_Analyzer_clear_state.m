function analyzer = Gfb_Analyzer_clear_state(analyzer)
% analyzer = Gfb_Analyzer_clear_state(analyzer)
%
% method for resetting the filter states to zeros
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006, Feb 2007

% filename : Gfb_Analyzer_clear_state.m

  for band = [1:length(analyzer.center_frequencies_hz)]
    analyzer.filters(1, band) = ...
	Gfb_Filter_clear_state(analyzer.filters(1, band));
  end
