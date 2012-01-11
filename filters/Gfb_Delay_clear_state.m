function delay = Gfb_Delay_clear_state(delay)
% delay = Gfb_Delay_clear_state(delay)
%
% Gfb_Delay_clear_state returns a copy of delay with cleared delay lines.
%
% PARAMETER:
% delay    A Gfb_Delay structure as returned by Gfb_Delay_new.  The returned
%          copy will have cleared delay lines
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_Delay_clear_state.m


delay.memory = zeros(size(delay.memory));

%OLDFORMAT
