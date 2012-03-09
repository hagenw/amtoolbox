function bp = amtbasepath;
%AMTBASEPATH  The base path of the AMT installation
%   Usage: bp = amtbasepath;
%
%   `amtbasepath` returns the top level directory in which the AMT
%   files are installed.
%
%   See also: amtstart
  
f=mfilename('fullpath');

bp = f(1:end-11);
