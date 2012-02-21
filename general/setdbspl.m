function inoutsig = setdbspl(inoutsig,lvl,varargin);
%SETDBSPL  Set level of signal in dB
%   Usage: outsig = setlevel(insig,lvl);
%          outsig = setlevel(insig,lvl,'ac');
%
%   `setdbspl(insig,lvl)` sets the SPL (sound pressure level) of the signal
%   insig to *lvl* dB, using the convention that a pure tone with an RMS value
%   of 1 corresponds to 100 dB SPL.
%
%   `setdbspl(lvl)` returns a scaling constant that will scale a signal
%   with RMS value of 1 to the correct level. 
%
%   If the input is a matrix, it is assumed that each column is a signal.
%
%   `setdbspl(insig,lvl,'ac')` does the same, but considers only the AC
%   component of the signal (i.e. the mean is removed).
%
%   References: moore2003introduction
  
%   See also: dbspl
  
%   Author: Peter L. Soendergaard, 2009

% ------ Checking of input parameters ---------
  
error(nargchk(1,3,nargin));

if ~isnumeric(inoutsig)
  error('%s: insig must be numeric.',upper(mfilename));
end;

% In the code below, "setdbspl" obtains the reference level from "dbspl"
% by calling "dbspl(1)", which will return only the offset measured in dB.

if nargin==1
  % Special mode, only the level has been given
  lvl=inoutsig;
  
  if ~isscalar(lvl) 
    error('%s: lvl must be a scalar.',upper(mfilename));
  end;

  inoutsig=gaindb(1,lvl-dbspl(1));
  return;
end;


if ~isnumeric(lvl) || ~isscalar(lvl) 
  error('%s: lvl must be a scalar.',upper(mfilename));
end;

if (nargin<3) || (~ischar(options))
  options='';
end;


% ------ Computation --------------------------

if isvector(inoutsig)
  inoutsig = gaindb(inoutsig/rms(inoutsig,varargin{:}),lvl-dbspl(1));
else
  % If we have a matrix, set the level for every column.
  for ii=1:size(inoutsig,2);
    inoutsig(:,ii) = gaindb(inoutsig(:,ii)/rms(inoutsig(:,ii),varargin{:}),lvl-dbspl(1));
  end;
end;