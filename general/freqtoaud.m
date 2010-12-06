function aud = freqtoaud(freq,varargin);
%FREQTOAUD  Converts frequencies (Hz) to auditory scale units.
%   Usage: aud = freqtoaud(scale,freq);
%
%   FREQTOAUD(freq,scale) converts values on frequency scale (measured in Hz) to
%   values on the selected auditory scale. The value of the parameter
%   scale determines the auditory scale:
%
%-    'mel'  -  The mel scale is a perceptual scale of pitches judged by
%               listeners to be equal in distance from one another. The
%               reference point between this scale and normal frequency
%               measurement is defined by equating a 1000 Hz tone, 40 dB above
%               the listener's threshold, with a pitch of 1000 mels.
%               The mel-scale is defined in Stevens et. al (1937).
%
%-    'bark'  - The bark-scale is originally defined in Zwicker (1961). A
%               distance of 1 on the bark scale is known as a critical
%               band. The implementation provided in this function is
%               described in Traunmuller (1990).
%
%-    'erb'   - A distance of 1 erb is equal to the equivalent rectangular
%               bandwidth of the auditory filters at that point on the
%               frequency scale. The scale is normalized such that 0 erbs
%               corresponds to 0 Hz. The width of the auditory filters were
%               determined by a notched-noise experiment. The erb scale is
%               defined in Glasberg and Moore (1990).
%
%-    'erb83' - This is the original defintion of the erb scale given in
%               Moore. et al. (1983).
%
%   If no flag is given, the function will print the list of valid flags.
%
%   See also: freqtoaud, audspace, audfiltbw
%
%R  stevens1937smp zwicker1961saf glasberg1990daf traunmuller1990aet moore1983sfc
  
%   AUTHOR: Peter L. Soendergaard

%% ------ Checking of input parameters ---------

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(freq) ||  all(freq(:)<0)
  error('%s: freq must be a non-negative number.',upper(mfilename));
end;

definput.flags.scale={'missingflag','mel','bark','erb','erb83'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.scale{2:end-2}),...
             sprintf('%s or %s',definput.flags.scale{end-1},...
                     definput.flags.scale{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% ------ Computation --------------------------


if flags.do_mel
  aud = 1127.01048*log(1+freq/700);
end;

if flags.do_erb
  aud = 9.265*log(1+freq/228.8455);
end;

if flags.do_bark
  % The bark scale seems to have several different approximations available.
  
  % This one was found through http://www.ling.su.se/STAFF/hartmut/bark.htm
  aud = (26.81./(1+1960./freq))-0.53;
  
  % The one below was found on Wikipedia.
  %aud = 13*atan(0.00076*freq)+3.5*atan((freq/7500).^2);
end;

if flags.do_erb83
  aud = 11.17*log((freq+312)./(freq+14675))+43.0;
end;
