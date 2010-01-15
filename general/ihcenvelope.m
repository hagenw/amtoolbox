function inoutsig = ihcenvelope(inoutsig,fs,methodname)
%IHCENVELOPE   Inner hair cell envelope extration
%   Usage:  outsig=ihcenvelope(insig,fs,methodname);
%
%   IHCENVELOPE(insig,fs,methodname) extract the envelope of an input signal
%   insig sampled with a sampling frequency of fs Hz. The envelope
%   extraction is performed by half-wave rectification followed by low pass
%   filtering. This is a common model of the signal transduction of the
%   inner hair cells.
%
%   The parameter methodname describes the kind of low pass filtering to
%   use. The name refers to a set of papers where in this particular
%   method has been utilized or studied. The options are
%
%-    'dau'     - Use a 2nd order butterworth filter with a cut-off
%                 frequency of 1000 Hz. This method has been used in all
%                 models deriving from the original 1996 model by 
%                 Dau et. al. These models are mostly monaural in nature.
%
%-    'hilbert' - Use the Hilbert envelope instead of the half-wave
%                 rectification and low pass filtering. This is not a
%                 releastic model of the inner hair envelope extraction
%                 process, but the option is included for completeness.  
%
%R  dau1996qmeI
  
% FIXME: Which paper did this idea originally appear in?

% ------ Checking of input parameters --------------------------------

error(nargchk(3,3,nargin));

if ~isnumeric(inoutsig)
  error('%s: The input signal must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

if ~ischar(methodname) 
  error('%s: methodname must be a string.',upper(mfilename));
end;

% ------ Computation -------------------------------------------------

if ~strcmp(lower(methodname),'hilbert')
  % Half-wave rectification
  inoutsig = max( inoutsig, 0 );
end;

switch(lower(methodname))
 case 'dau'
   cutofffreq=1000;
   [b, a] = butter(2, cutofffreq*2/fs);
   inoutsig = filter(b,a, inoutsig);
 case 'hilbert'
  inoutsig = abs(hilbert(inoutsig));
 otherwise
  error('%s: Unknown method name: %s.',upper(mfilename),methodname);
end;