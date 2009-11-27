function [y,n] = erbspacebw(varargin)
%ERBSPACE  Erbscale points specified by bandwidth
%   Usage: y=erbspacebw(flow,fhigh,bw,hitme);
%          y=erbspacebw(flow,fhigh,bw);
%          y=erbspacebw(flow,fhigh);
%
%   This is a wrapper around audspacebw that selects the erb-scale. Please
%   see the help on AUDSPACEBW for more information.
%M
%   See also: audspacebw, freqtoaud
  
%   AUTHOR : Peter L. Soendergaard
  
[y,n] = audspacebw('erb',varargin{:});