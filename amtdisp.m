function amtdisp(msg,flag)
%AMTDISP AMT-specific overload of the function 'disp'
%   Usage: amtdisp(X);
%     amtdisp(X,'progress');
%     amtdisp(X,'volatile');
%
%   `amtdisp(X);` can be used to show message X in the command
%   window. The output of `amtdisp` depends on the start-up
%   configuration of the AMT. 
%     
%   When the AMT is started in the `verbose` mode (default mode), `amtdisp` 
%   will always display. When the AMT is started in the `documentation` mode, 
%   `amtdisp` will display only if no flag is provided. When the AMT is 
%   started in the `silent` mode, `amtdisp` will never display. See
%   `amtstart` for further explanation on the start-up configurations. 
%
%   `amtdisp(X,'progress');` can be used as progress indicator. It will be
%   shown during the normal operation but supressed when used to create the
%   documentation.
%
%   `amtdisp(X,'volatile');` can be used as volatile progress indicator.
%   Any subsequent call of the `amtdisp` will delete the previous `volatile`
%   message. This way a changing progress can be clearly shown even in loops. 
% 
%   See also: amtstart

  
%   Author: Piotr Majdak, 2016
%   last change: 17.5.2017, volatile added

persistent CachedMsg;

if exist('flag','var')
  if ~strcmp(flag,'progress') && ~strcmp(flag,'volatile'),
    error(['Unsupported flag ' flag]);
  end
else
  flag='';
end

flags=amtflags;

if flags.do_verbose
  if strcmp(flag,'volatile')
    if ~isempty(CachedMsg),
      reversemsg = repmat(sprintf('\b'), 1, length(CachedMsg));
      fprintf(reversemsg);
    end
    fprintf(strrep(msg,'\','\\'));
    CachedMsg=msg;
  else
    if ~isempty(CachedMsg),
      reversemsg = repmat(sprintf('\b'), 1, length(CachedMsg));
      fprintf(reversemsg);
      CachedMsg=[];
    end
    disp(msg);
  end
end

if flags.do_documentation
  if ~strcmp(flag,'progress'), 
    disp(msg); 
  end
end

if flags.do_silent
  % do nothing
end