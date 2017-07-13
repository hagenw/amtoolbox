function amt_disp(msg,flag)
%amt_disp AMT-specific overload of the function 'disp'
%   Usage: amt_disp(X);
%     amt_disp(X,'progress');
%     amt_disp(X,'volatile');
%
%   `amt_disp(X);` can be used to show message X in the command
%   window. The output of `amt_disp` depends on the start-up
%   configuration of the AMT. 
%     
%   When the AMT is started in the `verbose` mode (default mode), `amt_disp` 
%   will always display. When the AMT is started in the `documentation` mode, 
%   `amt_disp` will display only if no flag is provided. When the AMT is 
%   started in the `silent` mode, `amt_disp` will never display. See
%   `amt_start` for further explanation on the start-up configurations. 
%
%   `amt_disp(X,'progress');` can be used as progress indicator. It will be
%   shown during the normal operation but supressed when used to create the
%   documentation.
%
%   `amt_disp(X,'volatile');` can be used as volatile progress indicator.
%   Any subsequent call of the `amt_disp` will delete the previous `volatile`
%   message. This way a changing progress can be clearly shown even in loops. 
% 
%   See also: amt_start

  
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

flags=amt_flags;

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