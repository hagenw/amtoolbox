function synthesizer = gfb_synthesizer_clear_state(synthesizer)
% synthesizer = gfb_synthesizer_clear_state(synthesizer)
%
% gfb_synthesizer_clear_state returns a copy of the synthesizer argument with
% a cleared state. (The delay lines will be cleared)
%
%  Input parameters:
% synthesizer  A gfb_Synthesizer structure as returned by gfb_synthesizer_new
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : gfb_synthesizer_clear_state.m

warning('Warning: GFB_SYNTHESIZER_CLEAR_STATE will be removed in a future release. Use HOHMANN2002CLEARSTATE instead. ');


synthesizer.delay = gfb_delay_clear_state(synthesizer.delay);


%OLDFORMAT
