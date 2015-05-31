function doit=amtredofile(filename,varargin)
%AMTREDOFILE  Determine if file should be redone
%   Usage: doit=amtredofile(filename,mode);
%
%   `amtredofile(file,mode)` returns 1 if the file should be redone, based
%   on the setting described in `mode`, otherwise is returns 0.
%
%   `mode` may one of the following flags:
%
%     'autorefresh'  Re-calculate the file if it does not exist. Return 1 if the
%                    file exist, otherwise 0. This is the default
%
%     'refresh'      Always recalculate the file.
%                  
%     'cached'       Always use the cached version. Throws an error if the
%                    file does not exist.

warning('amtredofile is deprecated and will be removed in future releases.\n For caching your data use AMTCACHE instead.');
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'amtredofile'};
[flags,kv]=ltfatarghelper({},definput,varargin);

doit=1;

if flags.do_refresh
  return;
end;

if exist(filename,'file')
  doit=0;  
else
  if flags.do_cached
    f=dbstack;  
    callfun=f(2).name;
    error('%s: A cached version of %s was requested, but it does not exist.',upper(callfun),filename);
  end;
end;