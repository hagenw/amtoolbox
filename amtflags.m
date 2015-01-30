function [flags,kv]=amtflags(varargin)
%AMTFLAGS Returns the start-up flags of the AMT
  
%   AUTHOR : Piotr Majdak 

persistent AMT;

  definput.flags.disp={'verbose','documentation','silent'};
%   definput.flags.cache={'download','redo','redowhenmissing'};

if ~isempty(varargin),
    % create persistent variable with flags
  if iscell(varargin{1}), varargin=varargin{1}; end
  [flags,kv]=ltfatarghelper({},definput,varargin);
  AMT.flags=flags;
  AMT.kv=kv;
elseif isempty(AMT)
  [flags,kv]=ltfatarghelper({},definput,[]);
  AMT.flags=flags;
  AMT.kv=kv;  
end
flags=AMT.flags;
kv=AMT.kv;
