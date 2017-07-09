function [hartestim,fs]  = data_harte2009()
%DATA_HARTE2009 Tone burst stimuli from Harte et al. (2009)
%   Usage: [hartestim,fs] = data_harte2009;
%
%   `[hartestim,fs]=data_harte2009` returns the tone burst stimuli from
%   Harte et al. (2009) and the sampling frequency, $fs=48000$.
%
%   References: harte2009comparison

% TODO: explain stimuli in description;

hartestim = amt_load('harte2009','stim.mat');
  
fs = 48e3;
  