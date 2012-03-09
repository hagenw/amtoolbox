function [ur,fs]  = data_roenne2012(varargin)
%DATA_ROENNE2012 Unitary response
%   Usage: ur = data_roenne2012()
%
%   `[ur,fs]=data_roenne2012` returns the unitary response from Roenne
%   (2012) and its sampling frequency, $fs=30000$.
%
%   References: roenne2012modeling
  
s=load([amtbasepath,'humandata',filesep,'ur']);
ur=s.ur;

fs=30000;


