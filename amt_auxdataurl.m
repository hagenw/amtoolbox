function auxURL=amt_auxdataurl(newURL)
%amt_auxdataurl URL of the auxiliary data
%   Usage: url=amt_auxdataurl
%          amt_auxdataurl(newurl)
%
%   `url=amt_auxdataurl` returns the URL of the web address containing
%   auxiliary data.
% 
%   `amt_auxdataurl(newurl)` sets the URL of the web address for further calls
%   of `amt_auxdataurl`.
%
%   See also: amt_auxdatapath amt_load

persistent AuxDataURL;

if exist('newURL','var')
  AuxDataURL=newURL;
elseif isempty(AuxDataURL)
  AuxDataURL=['http://www.sofacoustics.org/data/amt-' amt_help('version') '/auxdata'];
end
auxURL=AuxDataURL;
