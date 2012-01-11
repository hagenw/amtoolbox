function data = data_zwicker1961()
%DATA_ZWICKER1961  Data for the Bark scale
%   Usage: data = data_zwicker1961;
%
%   DATA_ZWICKER1961 returns the data points defining the notion of
%   critical bands. The output data consists of the frequency limit of
%   each band in Hz.
%
%   To get the bandwidth of each channel, simply use
%
%C     bw = diff(data_zwicker1961);
%
%   The first entry has been modified from the original paper: It was
%   originally 20 Hz, but in Zwicker and Fastl 1999 this has been changed
%   to 0 Hz.
%  
%   References:zwicker1961saf zwicker1999psychoacoustics
  
data = [0,100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,...
        2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500].';

%OLDFORMAT
