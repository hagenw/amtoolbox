function amtdisp(in,flag)
%AMTDISP AMT-specific overload of the function 'disp'
%   Usage: amtamtdisp(X);
%     amtdisp(X,'progress');
%
%   `amtdisp(X);` can be used for displaying information in the command
%   window in AMT functions. The output of amtdisp depends on the start-up
%   configuration of the AMT. 
%     
%   When the AMT is started in the 'verbose' mode, amtdisp will always
%   display. 
%
%   When the AMT is started in the 'documentation' mode, amtdisp will
%   display unless supressed by the flag 'progress' is provided. Thus, 
%   `amtdisp(in,'progress');` can be used as progress indicator when 
%   used in interactive way with the user but supress the progress in the
%   documentation.
%
%   When the AMT is started in the 'silent' mode, amtdisp will never
%   display. 
  
%   Author: Piotr Majdak, 2014

if exist('flag','var')
  if ~strcmp(flag,'progress'),
    error(['Unsupported flag ' flag]);
  end
else
  flag='';
end

flags=amtflags;

if flags.do_verbose
  disp(in);
end

if flags.do_documentation
  if ~strcmp(flag,'progress'), disp(in); end
end

if flags.do_silent
  % do nothing
end