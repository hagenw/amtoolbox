function varargout=amtcache(cmd,name,varargin)
%AMTCACHE  Cache variables for later or retrieves variables from cache
%   Usage: var = amtcache('get',package,flags);
%     amtcache('set',package,variables);
%   
%   `amtcache` supports the following commands:
%
%     'get'      gets the content of a package from the cache. 
%                `variables = amtcache('get',package)` reads a `package` from the cache
%                and outputs its content in `var`. `package` must a be a string
%                identifying the package of variables. If the package contains multiple
%                variables, `variables` can be a list of variables like 
%                `[var1, var2, ... , varN] = ...`. The order of returned variables is the
%                same as that used for saving in cache.
%                `... = amtcache('get',package,flags)` allows to control the 
%                behaviour of accessing the cache. `flags` can be:
%                  'normal': package will be recalculated when locally not available. 
%                  'redo': enforce the recalculation of the package. [..] = amtcache('get', [..]) outputs empty variables always. 
%                  'cached': enforce using cached package. If the cached package is locally not available, it will be downloaded from the internet. If it is remotely not available, warning will be thrown and the package will be recalculated. Note that this method may by-pass the actual processing and thus does not test the corresponding functionality. It is, however, very convenient for fast access of results like plotting figures. On the internet, the cached packages are available only for the models from release version of the AMToolbox. 
%
%     'set'      stores variables as a package in the cache. 
%                `amtcache('set',package, variables)` saves variables in the cache using
%                the name `package`. `variables` can be a list of variables separated by
%                comma.
%                
%     'getURL'   outputs the URL of the cache in the internet. 
%
%     'setURL'   sets the URL of the internet cache to a new URL. 
%
%     'clearAll' clears the cache directory. An interactive confirmation is
%                required.
%
%
%   This is an example of using the cache in a function:
%
%     definput.import={'amtcache'};
%     [flags,~]  = ltfatarghelper({},definput,varargin);
%
%     [x,y,z] = amtcache('get', 'xyz', flags.cachemode);
% 
%     if isempty(x)
%         calculate your variables x,y,z here
%         amtcache('set','xyz',x,y,z);
%     end
%     use your variables x,y,z here
%
%   See also: data_ziegelwanger2013



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