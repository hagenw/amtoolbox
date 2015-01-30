function auxURL=amtauxdataurl(newURL)
%amtauxdataurl URL of the auxiliary data
%   Usage: url=amtauxdataurl
%          amtauxdataurl(newurl)
%
%   `url=amtauxdataurl` returns the URL of the web address containing
%   auxiliary data.
% 
%   `amtauxdataurl(newurl)` sets the URL of the web address for further calls
%   of `amtauxdataurl`.
%
%   See also: amtauxdatapath amtload

persistent CachedURL;

if exist('newURL','var')
  CachedURL=newURL;
elseif isempty(CachedURL)
  CachedURL=['http://www.sofacoustics.org/data/amt-' amthelp('version') '/auxdata'];
end
auxURL=CachedURL;
