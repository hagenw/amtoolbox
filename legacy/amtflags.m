function [flags,kv]=amtflags(varargin)

warning('Warning: AMTFLAGS will be removed in a future release. Use AMT_FLAGS instead. ');  

% Based on amt_flags from 9.7.2017

persistent AMT;
  
if ~isempty(varargin),
    % create persistent variable with flags
  if iscell(varargin{1}), varargin=varargin{1}; end
  definput.import={'amt_cache','amt_disp'};
  [flags,kv]=ltfatarghelper({},definput,varargin);
  AMT.flags=flags;
  AMT.kv=kv;
elseif isempty(AMT)
  definput.import={'amt_cache','amt_disp'};
  [flags,kv]=ltfatarghelper({},definput,[]);
  AMT.flags=flags;
  AMT.kv=kv;  
end
flags=AMT.flags;
kv=AMT.kv;
