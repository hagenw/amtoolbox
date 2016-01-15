function [output, synthesizer] = gfb_synthesizer_process(synthesizer, input)
% [output, synthesizer] = gfb_synthesizer_process(synthesizer, input)
%
% The synthesizer will resynthesize the given input.
%
%  Input parameters:
% synthesizer  A synthesizer structure as created by gfb_synthesizer_new. A
%              copy of the synthesizer object with an updated internal state
%              is returned in the second return parameter
% input        A matrix containing the (possibly processed) complex output of
%              the analyzer corresponding to this synthesizer.  The number of
%              rows in input must match the number of filter bands
% output       The synthesized output signal
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : gfb_synthesizer_process.m

warning('Warning: GFB_SYNTHESIZER_PROCESS will be removed in a future release. Use HOHMANN2002PROCESS instead. ');

[output, synthesizer.delay] = gfb_delay_process(synthesizer.delay, input);
[output, synthesizer.mixer] = gfb_mixer_process(synthesizer.mixer, output);

%OLDFORMAT
