function varargout=amtcache(cmd,name,varargin)
%AMTCACHE  Cache variables for later or retrieves variables from cache
%   Usage: amtcache(..);
%
%  


%   Author: Piotr Majdak, 2015

persistent CacheURL;
if isempty(CacheURL)
  CacheURL=['http://www.sofacoustics.org/data/amt-' amthelp('version') '/cache'];
end

switch cmd
  case 'set'
    f=dbstack('-completenames');
    fn=f(2).file;
    token=urlencode(strrep(fn(length(amtbasepath)+1:end),'\','/'));
    tokenpath=fullfile(amtbasepath,'cache',token);
    tokenfn=fullfile(tokenpath,[name '.mat']);
    
    if ~exist(tokenpath,'dir'); mkdir(tokenpath); end
    for ii=3:nargin
      cache(ii-2).name=inputname(ii);
      cache(ii-2).value=varargin{ii-2};
    end
    save(tokenfn,'cache','-v6');
    varargout{1}=tokenfn;
  case 'get'
    if nargin<3, varargin{1}='normal'; end    
    switch varargin{1}
      case 'redo' % force recalculations in any case
        for ii=1:nargout, varargout{ii}=[]; end

      case 'cached' % force using cached in any case
        f=dbstack('-completenames');
        fn=f(2).file;
        token=urlencode(strrep(fn(length(amtbasepath)+1:end),'\','/'));
        tokenpath=fullfile(amtbasepath,'cache',token);
        tokenfn=fullfile(tokenpath,[name '.mat']);

        if ~exist(tokenfn,'file'),
          % TODO: implement downloading from the internet
        end
        load(tokenfn);
        for ii=1:nargout
          varargout{ii}=s(ii).value;
        end
        
      case {'normal','global'} % normal mode: redo if locally not available
        f=dbstack('-completenames');
        fn=f(2).file;
        token=urlencode(strrep(fn(length(amtbasepath)+1:end),'\','/'));
        tokenpath=fullfile(amtbasepath,'cache',token);
        tokenfn=fullfile(tokenpath,[name '.mat']);

        if exist(tokenfn,'file'),
          load(tokenfn);
          for ii=1:nargout
            varargout{ii}=cache(ii).value;
          end
        else
          for ii=1:nargout, varargout{ii}=[]; end
        end            
    end
  case 'setURL'
    CacheURL=name;
  case 'getURL'
    varargout{1}=CacheURL;
  case 'clearAll'
    cachepath=fullfile(amtbasepath,'cache');
    if strcmp(input(['clearAll clears ' strrep(cachepath,'\','\\') '. Type YES for confirmation: '],'s'),'YES'), 
      amtdisp(['Clearing ' cachepath ' ...']);
      rmdir(cachepath, 's');       
    end
  otherwise
    error('Unsupported command');
end