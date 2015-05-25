function varargout=amtcache(cmd,name,varargin)
%AMTCACHE  Cache variables for later or retrieves variables from cache
%   Usage: amtcache(..);
%
%  


%   Author: Piotr Majdak, 2015


f=dbstack('-completenames');
fn=f(2).file;

token=urlencode(strrep(fn(length(amtbasepath)+1:end),'\','/'));
tokenpath=fullfile(amtbasepath,'cache',token);
tokenfn=fullfile(tokenpath,[name '.mat']);

switch cmd
  case 'set'
    if ~exist(tokenpath,'dir'); mkdir(tokenpath); end
    for ii=3:nargin
      s(ii-2).name=inputname(ii);
      s(ii-2).value=varargin{ii-2};
    end
    save(tokenfn,'s','-v6');
    varargout{1}=tokenfn;
  case 'get'
    if nargin<3, varargin{1}='cached'; end
    switch varargin{1}
      case 'refresh'
        for ii=1:nargout, varargout{ii}=[]; end
      case {'autorefresh','cached'}
        if exist(tokenfn,'file'),
          load(tokenfn);
          for ii=1:nargout
            varargout{ii}=s(ii).value;
          end
        else
          for ii=1:nargout, varargout{ii}=[]; end
        end    
    end
  case 'delete'
    varargout{1}=[];
  otherwise
    varargout{1}=[];
end