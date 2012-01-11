function filter = Gfb_Filter_clear_state(filter)
% filter = Gfb_Filter_clear_state(filter)
% 
% returns a copy of the filter, with the filter state cleared.
%
% PARAMETER:
% filter  a Gfb_Filter structure as returned by Gfb_Filter_new.  A copy
%         of the filter is returned, with the filter state cleared
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_Filter_clear_state.m


filter.state = zeros(1, filter.gamma_order);


%OLDFORMAT
