function [output, mixer] = gfb_mixer_process(mixer, input)
% [output, mixer] = gfb_mixer_process(mixer, input)
%
% gfb_mixer_process computes a weighted sum of the different bands present
% in input.
%
%  Input parameters:
% mixer   A gfb_mixer structure as returned by gfb_mixer_new.  The mixer
%         contains the gain factors for the weighted sum.  A copy of mixer
%         will be returned in the second return parameter
% input   an NxM matrix, where N equals the number of gain factors (bands)
%         of the mixer
% output  an 1xM vector containing the weighted sums of each comlumn
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : gfb_mixer_process.m

warning('Warning: GFB_MIXER_PROCESS will be removed in a future release. Use HOHMANN2002PROCESS instead. ');

output = mixer.gains * input;


%OLDFORMAT
