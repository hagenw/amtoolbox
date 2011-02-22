function inoutsig = setdbspl(inoutsig,lvl,options);
%SETDBSPL  Set level of signal in dB
%   Usage: outsig = setlevel(insig,lvl);
%          outsig = setlevel(insig,lvl,'ac');
%
%   SETDBSPL(insig,lvl) sets the SPL (sound pressure level) of the signal
%   insig to lvl dB, using the convention that a pure tone with an RMS value
%   of 1 corresponds to 100 dB SPL.
%
%   If the input is a matrix, it is assumed that each column is a signal.
%
%   SETDBSPL(insig,lvl,'ac') does the same, but considers only the AC
%   component of the signal (i.e. the mean is removed).
%
%R  moore2003introduction
  
%   See also: dbspl
  
%   Author: Peter L. Soendergaard, 2009

% ------ Checking of input parameters ---------
  
error(nargchk(2,3,nargin));

if ~isnumeric(inoutsig)
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(lvl) || ~isscalar(lvl) 
  error('%s: lvl must be a scalar.',upper(mfilename));
end;

if (nargin<3) || (~ischar(options))
  options='';
end;

% ------ Computation --------------------------

if isvector(inoutsig)
  inoutsig = gaindb(inoutsig/rms(inoutsig,options),lvl-100);
else
	% If we have a matrix, set the level for every column.
  for ii=1:size(inoutsig,2);
    inoutsig(:,ii) = gaindb(inoutsig(:,ii)/rms(inoutsig(:,ii),options),lvl-100);
  end;
end;
