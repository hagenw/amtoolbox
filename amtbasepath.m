function bp = amtbasepath;
%AMTBASEPATH  The base path of the AMT installation
%   Usage: bp = amtbasepath;
%
%   AMTBASEPATH returns the top level directory in which the AMT
%   files are installed.
%
%   See also: amtstart
  
global AMT_CONF

bp = AMT_CONF.basepath;