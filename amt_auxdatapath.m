function auxPath=amt_auxdatapath(newPath)
%amt_auxdatapath Local path to the auxiliary data
%   Usage: auxpath=amt_auxdatapath
%          amt_auxdatapath(newpath)
%
%   `auxPath=amt_auxdatapath` returns the path of the directory containing
%   auxiliary data.
%
%   Default path to the auxiliary data is the `amt_basepath/auxdata`.
% 
%   `amt_auxdatapath(newpath)` sets the path of the directory for further calls
%   of `amt_auxdatapath`.
%
%   See also: amt_auxdataurl amt_load amt_basepath

% AUTHOR: Piotr Majdak, 2015


persistent CachedPath;

if exist('newPath','var')
  CachedPath=newPath;
elseif isempty(CachedPath)
  CachedPath=fullfile(amt_basepath, 'auxdata');
end
auxPath=CachedPath;

  