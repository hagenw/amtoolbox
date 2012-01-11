function [output, synthesizer] = Gfb_Synthesizer_process(synthesizer, input)
% [output, synthesizer] = Gfb_Synthesizer_process(synthesizer, input)
%
% The synthesizer will resynthesize the given input.
%
% PARAMETERS:
% synthesizer  A synthesizer structure as created by Gfb_Synthesizer_new. A
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

% filename : Gfb_Synthesizer_process.m


[output, synthesizer.delay] = Gfb_Delay_process(synthesizer.delay, input);
[output, synthesizer.mixer] = Gfb_Mixer_process(synthesizer.mixer, output);

%OLDFORMAT
