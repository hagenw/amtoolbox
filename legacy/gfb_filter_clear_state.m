function filter = gfb_filter_clear_state(filter)
%GFB_FILTER_CLEAR_STATE  Clear filters
%   Usage: filter = gfb_filter_clear_state(filter)
%
%   Input parameters:
%     filter : a `gfb_filter` structure as returned by
%              |gfb_filter_new|.
%
%   `gfb_filter_clear_state(filter)` returns a copy of the filter, with the
%   filter state cleared.
%
%   See also: gfb_filter_new

% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

warning('Warning: GFB_FILTER_CLEAR_STATE will be removed in a future release. Use HOHMANN2002CLEARSTATE instead. ');

filter.state = zeros(1, filter.gamma_order);
