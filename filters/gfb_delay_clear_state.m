function delay = gfb_delay_clear_state(delay)
% delay = gfb_delay_clear_state(delay)
%
% gfb_delay_clear_state returns a copy of delay with cleared delay lines.
%
%  Input parameters::
% delay    A gfb_Delay structure as returned by gfb_delay_new.  The returned
%          copy will have cleared delay lines
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : gfb_delay_clear_state.m


delay.memory = zeros(size(delay.memory));

%OLDFORMAT
