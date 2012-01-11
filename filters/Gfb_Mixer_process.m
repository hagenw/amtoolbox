function [output, mixer] = Gfb_Mixer_process(mixer, input)
% [output, mixer] = Gfb_Mixer_process(mixer, input)
%
% Gfb_Mixer_process computes a weighted sum of the different bands present
% in input.
%
% PARAMETERS:
% mixer   A Gfb_Mixer structure as returned by Gfb_Mixer_new.  The mixer
%         contains the gain factors for the weighted sum.  A copy of mixer
%         will be returned in the second return parameter
% input   an NxM matrix, where N equals the number of gain factors (bands)
%         of the mixer
% output  an 1xM vector containing the weighted sums of each comlumn
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_Mixer_process.m


output = mixer.gains * input;


%OLDFORMAT
