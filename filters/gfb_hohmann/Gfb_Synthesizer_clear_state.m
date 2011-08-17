function synthesizer = Gfb_Synthesizer_clear_state(synthesizer)
% synthesizer = Gfb_Synthesizer_clear_state(synthesizer)
%
% Gfb_Synthesizer_clear_state returns a copy of the synthesizer argument with
% a cleared state. (The delay lines will be cleared)
%
% PARAMETER
% synthesizer  A Gfb_Synthesizer structure as returned by Gfb_Synthesizer_new
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_Synthesizer_clear_state.m


synthesizer.delay = Gfb_Delay_clear_state(synthesizer.delay);

