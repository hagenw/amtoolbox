function [flags,kv]=amt_flags(varargin)
%amt_flags Returns the start-up flags of the AMT
  
%   AUTHOR : Piotr Majdak 

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
