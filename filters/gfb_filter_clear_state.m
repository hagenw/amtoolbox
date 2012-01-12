function filter = gfb_filter_clear_state(filter)
% filter = gfb_filter_clear_state(filter)
% 
% returns a copy of the filter, with the filter state cleared.
%
%  Input parameters::
% filter  a gfb_Filter structure as returned by gfb_filter_new.  A copy
%         of the filter is returned, with the filter state cleared
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : gfb_filter_clear_state.m


filter.state = zeros(1, filter.gamma_order);


%OLDFORMAT
