function [y,bw] = audspace(scale,flow,fhigh,n)
%AUDSPACE  Equidistantly spaced points on auditory scale
%   Usage: y=audspace(scale,flow,fhigh,n);
%
%   AUDSPACE(scale,flow,fhigh,n) computes a vector of length n containing values
%   equistantly scaled on the selected auditory scale between flow and fhigh. All
%   frequencies are specified in Hz.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter.
%  
%   [y,bw]=AUDSPACE(...) does the same but outputs the bandwidth between
%   each sample measured on the selected scale.
%  
%
%   See also: freqtoaud, audspacebw, audfiltbw
%
%R  glasberg1990daf
  
%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ---------
  
error(nargchk(4,4,nargin));
  
% Default parameters.

if ~isnumeric(flow) || ~isscalar(flow) || flow<0
  error('%s: flow must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) || fhigh<0
  error('%s: fhigh must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(n) || ~isscalar(n) || n<=0 || fix(n)~=n
  error('%s: n must be a positive, integer scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;


% ------ Computation --------------------------

audlimits = freqtoaud(scale,[flow,fhigh]);

y = audtofreq(scale,linspace(audlimits(1),audlimits(2),n));  

bw=(audlimits(2)-audlimits(1))/(n-1);
