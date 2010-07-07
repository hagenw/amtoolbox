function inoutsig = ihcenvelope(inoutsig,fs,varargin)
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
%                 by raising it to the power .2, combine the envelope
%                 with the original fine-structure, half-wave rectify it, 
%                 square it and low-pass filter it with a cut-off
%                 frequency of 425 Hz. This method is defined in
%                 Bernstein 1999. Note that this method includes both a
%                 compression and an expansion stage.
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
%                 process, but the option is included for
%                 completeness. The Hilbert envelope was first suggested
%                 for signal analysis in Gabor 1946.
%
%-    'lindemann' - Use a 1st order Butterworth filter with a cut-off frequency
%                 of 800 Hz. This method is defined in the paper
%                 Lindemann 1986a.
%
%-    'minlvl' - Set all values in the output equal to minlvl. The ensures that
%                the output is non-negative and that further processing
%                is not affected by unnaturally small values. The default
%                value of [] means to not do this.
%
%R  bernstein1999normalized breebaart2001binaural gabor1946 lindemann1986a dau1996qmeI
  
% FIXME: Which paper did this idea originally appear in?

% ------ Checking of input parameters --------------------------------

if nargin<2
  error('Too few input parameters.');
end;


if ~isnumeric(inoutsig)
  error('%s: The input signal must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.flags.model={'nodefault','bernstein','breebart','dau','hilbert', ...
                    'lindemann'};

definput.flags.keyvals.minlvl=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% ------ Computation -------------------------------------------------

if flags.do_nodefault
  error(['%s: you must supply a flag to designate the IHC model to ' ...
         'use.'],upper(mfilename));
end;

if flags.do_bernstein
  % The computational trick mentioned in the Bernstein paper is used
  % here: Instead of raising the envelope to power .23 and combine with its
  % TFS, we raise it to power -.77, and combine with the original
  % signal. In this way we avoid computing the fine structure.
  inoutsig=max(abs(hilbert(inoutsig)).^(-.77).*inoutsig,0).^2;
  cutofffreq=425;
  [b, a] = butter(2, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
end;

if flags.do_breebart
  inoutsig = max( inoutsig, 0 );
  cutofffreq=2000;
  [b, a] = butter(1, cutofffreq*2/fs);
  for ii=1:5
    inoutsig = filter(b,a, inoutsig);
  end;
end;

if flags.do_dau
  inoutsig = max( inoutsig, 0 );
  cutofffreq=1000;
  [b, a] = butter(2, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
end;

if flags.do_hilbert
  inoutsig = abs(hilbert(inoutsig));
end;

if flags.do_lindemann
  inoutsig = max( inoutsig, 0 );
  cutofffreq=800;
  [b, a] = butter(1, cutofffreq*2/fs);
  inoutsig = filter(b,a, inoutsig);
end;

if ~isempty(keyvals.minlvl)
  inoutsig = max( inoutsig, keyvals.minlvl );
end;

