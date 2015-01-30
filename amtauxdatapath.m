function auxPath=amtauxdatapath(newPath)
%amtauxdatapath Path of the auxiliary data
%   Usage: auxpath=amtauxdatapath
%          amtauxdatapath(newpath)
%
%   `auxPath=amtauxdatapath` returns the path of the directory containing
%   auxiliary data.
%
%   Default path to the auxiliary data is the `amtbasepath/auxdata`.
% 
%   `amtauxdatapath(newpath)` sets the path of the directory for further calls
%   of `amtauxdatapath`.
%
%   See also: amtauxdataurl amtload amtbasepath

% AUTHOR: Piotr Majdak, 2015


persistent CachedPath;

if exist('newPath','var')
  CachedPath=newPath;
elseif isempty(CachedPath)
  CachedPath=fullfile(amtbasepath, 'auxdata');
end
auxPath=CachedPath;

  