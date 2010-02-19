function [optargs,varargout]  = amtarghelper(nargs,defaults,varargin)
%AMTARGHELPER  Parse arguments for AMT
%   Usage: [flags,varargout]  = amt_arghelper(nargs,defaults,varargin);
%
%   Input parameters:
%      nargs     : Number of required parameters
%      defaults  : Cell array with default values
%      varargin  : Supply the commandline of the calling function
%
%   Output parameters:
%      optargs   : Cell array of optional arguments
%      varargout : The required paramers properly initialized
%
%   [optargs,varargout]=AMTARGHELP(nargs,defaults,varargin) assist in
%   parsing input parameters for a function in AMT. Parameters come in
%   two categories: the first given parameters must be numers of some
%   kind (not string, and if they are absent a default value is
%   given. The second category are optional parameters that come in
%   key,value pairs.
%
%   The typical way of calling AMTARGHELP is as follows:
%  
%C     [optargs,n,betamul] = amtarghelper(2,{4,[]},varargin{:});
%
%   This will pass any (or a number) of the input arguments from the
%   calling function onto AMTARGHELPER. In this case, there are 2
%   required arguments (n and betamul), which will have default values 4
%   and [] if they are not present.
  
  total_args = numel(varargin);
    
  % Determine the position of the first optional argument.
  % If no optional argument is given, return nargs+1
  first_str_pos = 1;
  while first_str_pos<=total_args && ~isstr(varargin{first_str_pos}) 
    first_str_pos = first_str_pos +1;    
  end;
  
  % If more than nargs arguments are given, the first additional one must
  % be a string
  if (first_str_pos>nargs+1)
    error('Too many input arguments');
  end;

  n_first_args=min(nargs,first_str_pos-1);
  
  % Copy the given first arguments
  for ii=1:n_first_args
    varargout(ii)=varargin(ii);
  end;

  % Copy the default values for the rest
  for ii=n_first_args+1:nargs
    varargout(ii)=defaults(ii);
  end;
  
  % Copy the optional arguments
  optargs={};
  for ii=0:total_args-first_str_pos
    optargs(ii+1)=varargin(first_str_pos+ii);
  end;
  
  
  
  