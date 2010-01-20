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
%-    'bernstein' - Compute the Hilbert envelope, compress the envelope
%                 by raising it to the power .23, half-wave rectify it, 
%                 square it and low-pass filter it with a cut-off
%                 frequency of 425 Hz. This method is defined in Bernstein 1999.
%
%-    'breebart' - Use a 5th order filter with a cut-off frequency of 770
%                 Hz. This method is given in Breebart 2001. Page 94 in thesis.
%
%-    'dau'     - Use a 2nd order Butterworth filter with a cut-off
%                 frequency of 1000 Hz. This method has been used in all
%                 models deriving from the original 1996 model by 
%                 Dau et. al. These models are mostly monaural in nature.
%
%-    'hilbert' - Use the Hilbert envelope instead of the half-wave
%                 rectification and low pass filtering. This is not a
%                 releastic model of the inner hair envelope extraction
%                 process, but the option is included for completeness.
%
%-    'lindemann' - Use a 1st order Butterworth filter with a cut-off frequency
%                 of 800 Hz. This method is defined in the paper
%                 Lindemann 1986a.
%
%R  bernstein1999normalized breebaart2001binaural lindemann1986a dau1996qmeI
  
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
 case 'bernstein'
  inoutsig=max(abs(hilbert(inoutsig)).^0.23,0).^2;
  cutofffreq=425;
  [b, a] = butter(2, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
 case 'breebart'
  inoutsig = max( inoutsig, 0 );
  cutofffreq=2000;
  [b, a] = butter(1, cutofffreq*2/fs);
  for ii=1:5
    inoutsig = filter(b,a, inoutsig);
  end;
 case 'dau'
  inoutsig = max( inoutsig, 0 );
  cutofffreq=1000;
  [b, a] = butter(2, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
 case 'hilbert'
  inoutsig = abs(hilbert(inoutsig));
 case 'lindemann'
  inoutsig = max( inoutsig, 0 );
  cutofffreq=800;
  [b, a] = butter(1, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
 otherwise
  error('%s: Unknown method name: %s.',upper(mfilename),methodname);
end;