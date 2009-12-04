function s=greasy()
%GREASY  Load the 'greasy' test signal.
%   Usage:  s=greasy;
%
%   GREASY loads the 'greasy' signal. It is a recording of a woman
%   pronouncing the word "greasy".
%
%   The signal is 5880 samples long and recorded at 16 khz.
%
%   The signal has been scaled to not produce any clipping when
%   played. To get the original integer values use round(greasy*1780).
%
%   The signal was obtained from Wavelab:
%     http://www-stat.stanford.edu/~wavelab/
%M
%R  mazh93

%   AUTHOR : Peter Soendergaard
  
if nargin>0
  error('This function does not take input arguments.')
end;

s=load('-ascii',[amtbasepath,'signals',filesep,'greasy.asc']);


