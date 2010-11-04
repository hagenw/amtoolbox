function freq = audtofreq(scale,aud);
%AUDTOFREQ  Converts auditory units to frequency (Hz)
%   Usage: freq = audtofreq(aud);
%  
%   AUDTOFREQ(scale,aud) converts values on the selected auditory scale to
%   values on the frequency scale measured in Hz.
%
%   See the help on FREQTOAUD to get a list of the supported values of the
%   scale parameter. 
%
%   See also: freqtoaud, audspace, audfiltbw

%   AUTHOR: Peter L. Soendergaard

% ------ Checking of input parameters ---------

error(nargchk(2,2,nargin));

if ~isnumeric(aud) ||  all(aud(:)<0)
  error('%s: aud must be a non-negative number.',upper(mfilename));
end;

if ~ischar(scale)
  error('%s: the scale must be denoted by a character string.',upper(mfilename))
end;

% ------ Computation --------------------------
  
switch(lower(scale))
 case 'mel'
  freq = 700*(exp(aud/1127.01048)-1);
 case 'erb'
  freq = 228.8455*(exp(aud/9.265)-1);
 case 'bark'
  % This one was found through http://www.ling.su.se/STAFF/hartmut/bark.htm
  freq = 1960./(26.81./(aud+0.53)-1);
 case 'erb83'
  freq = 14363./(1-exp((aud-43.0)/11.7))-14675;
 otherwise
  error(['%s: unknown auditory scale: %s. Please see the help for a list ' ...
         'of supported scales.'],upper(mfilename),scale);
end;