function [y,n] = audspacebw(scale,flow,fhigh,bw,hitme)
%AUDSPACEBW  Auditory scale points specified by bandwidth
%   Usage: y=audspacebw(scale,flow,fhigh,bw,hitme);
%          y=audspacebw(scale,flow,fhigh,bw);
%          y=audspacebw(scale,flow,fhigh);
%          [y,n]=audspacebw(...);
%
%   AUDSPACEBW(scale,flow,fhigh,bw,hitme) computes a vector containing
%   values equistantly scaled between flow and fhigh on the selected
%   auditory scale. The distance between two consecutive values are given by
%   a value of bw the scale. One of the points is quaranteed to be the
%   frequency hitme. All frequencies are specified in Hz.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter.
%  
%   AUDSPACEBW(scale,flow,fhigh,bw) will do the same but center the points between
%   flow and fhigh without including a specific frequency (this can also
%   be specified by setting a negative value for hitme).
%
%   AUDSPACEBW(scale,flow,fhigh) will do the same assuming a spacing of 1
%   on the selected scale.
%
%   [y,n]=AUDSPACEBW( ... ) additionally returns the number of points n in
%   the output vector y.
%
%   See also: freqtoaud, audspace, audfiltbw
  
%   AUTHOR : Peter L. Soendergaard
  
% ------ Checking of input parameters ---------
  
error(nargchk(3,5,nargin));
  
% Default parameters
if nargin<5
  hitme=-1;
end;

if nargin<4
  bw=1;
end;

if ~isnumeric(flow) || ~isscalar(flow) || flow<0
  error('%s: flow must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(fhigh) || ~isscalar(fhigh) || fhigh<0
  error('%s: fhigh must be a non-negative scalar.',upper(mfilename));
end;

if ~isnumeric(bw) || ~isscalar(bw) || bw<=0 
  error('%s: bw must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(hitme) || ~isscalar(hitme) 
  error('%s: hitme must be a scalar.',upper(mfilename));
end;

if flow>fhigh
  error('%s: flow must be less than or equal to fhigh.',upper(mfilename));
end;

if nargin==4 && hitme>=0
  if (hitme<flow) || (hitme>fhigh)
    error('%s: hitme must be in the interval from flow to fhigh.',upper(mfilename));
  end;
end;

% ------ Computation --------------------------

if hitme<0
  % Convert the frequency limits to auds.
  audlimits = freqtoaud(scale,[flow,fhigh]);
  audrange  = audlimits(2)-audlimits(1);

  % Calculate number of points, excluding final point
  n         = floor(audrange/bw);

  % The remainder is calculated in order to center the points
  % correctly between flow and fhigh.
  remainder = audrange-n*bw;

  audpoints = audlimits(1)+(0:n)*bw+remainder/2;
  
  % Add the final point
  n=n+1;
  
else
  % Convert the frequency limits to auds.
  audlimits    = freqtoaud(scale,[flow,fhigh,hitme]);
  audrangelow  = audlimits(3)-audlimits(1);
  audrangehigh = audlimits(2)-audlimits(3);

  % Calculate number of points, exluding final point.
  nlow = floor(audrangelow/bw);
  nhigh = floor(audrangehigh/bw);
  
  audpoints=(-nlow:nhigh)*bw+audlimits(3);
  n=nlow+nhigh+1;
end;

y = audtofreq(scale,audpoints);
