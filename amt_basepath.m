function bp = amt_basepath;
%amt_basepath  The base path of the AMT installation
%   Usage: bp = amt_basepath;
%
%   `amt_basepath` returns the top level directory in which the AMT
%   files are installed.
%
%   See also: amt_start amt_auxdatapath
  
f=mfilename('fullpath');

bp = f(1:end-12);
