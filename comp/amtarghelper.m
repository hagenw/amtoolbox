function [flags,keyvals,varargout]  = amtarghelper(nposdep,defposdep,defnopos,arglist,callfun)
%AMTARGHELPER  Parse arguments for AMT
%   Usage: [flags,varargout]  = amtarghelper(nposdep,defposdep,defnopos,arglist,callfun);
%
%   Input parameters:
%      nposdep   : Number of position dependant parameters
%      defposdep : Default values for pos. dep. pars. (cell array)
%      defnopos  : Definitions of the pos. independent pars (cell array)
%      arglist   : Commandline of the calling function
%      callfun   : Name of calling function
%
%   Output parameters:
%      nopos     : Cell array containing parsed pos. independent pars.
%      varargout : The position dependant pars. properly initialized
%
%   [nopos,varargout]=AMTARGHELPER(nposdep,defposdep,defnopos,arglist,callfun) 
%   assist in parsing input parameters for a function in AMT. Parameters come 
%   in two categories:
%  
%      * Position dependant parameters. These must not be strings. These are
%      the first parameters passed to a function. If they are not present a
%      default value is given.
%
%      * Position independant parameters. These come in key,value pairs
%      at the end of the parameter list of the calling function.
%
%   The typical way of calling AMTARGHELPER is as follows:
%  
%C     [nopos,n,betamul] = amtarghelper(2,{4,[]},defnopos,arglist{:},...
%           upper(mfilename));
%
%   This will pass any (or a number) of the input arguments from the calling
%   function onto AMTARGHELPER. In this case, there are 2 position
%   dependant parameters (n and betamul), which will have default values 4
%   and [] if they are not present.
%   defnopos is a struct and can have the fields flags and keyvals.

  if isfield(defnopos,'flags')
    defflags=defnopos.flags;
  else
    defflags=struct;
  end;

  if isfield(defnopos,'keyvals')
    defkeyvals=defnopos.keyvals;
  else
    defkeyvals=struct;
  end;

  
  
  total_args = numel(arglist);
    
  % Determine the position of the first optional argument.
  % If no optional argument is given, return nposdep+1
  first_str_pos = 1;
  while first_str_pos<=total_args && ~ischar(arglist{first_str_pos}) 
    first_str_pos = first_str_pos +1;    
  end;
    
  % If more than nposdep arguments are given, the first additional one must
  % be a string
  if (first_str_pos>nposdep+1)
    error('%s: Too many input arguments',callfun);
  end;

  n_first_args=min(nposdep,first_str_pos-1);
  
  % Copy the given first arguments
  for ii=1:n_first_args
    varargout(ii)=arglist(ii);
  end;

  % Copy the default values for the rest
  for ii=n_first_args+1:nposdep
    varargout(ii)=defposdep(ii);
  end;
  
  % Initialize the position independent parameters.
  % and create reverse mapping of flag -> group
  flagnames=fieldnames(defflags);
  flags=struct;
  flagsreverse=struct;
  for ii=1:numel(flagnames)
    name=flagnames{ii};
    flaggroup=defflags.(name);
    flags.(name)=flaggroup{1};
    for jj=1:numel(flaggroup)
      flagsreverse.(flaggroup{jj})=name;
      flags.(['do_',flaggroup{jj}])=0;
    end;
    flags.(['do_',flaggroup{1}])=1;
  end;
  
  keyvals=defkeyvals;
      
  ii=first_str_pos;
  while ii<=total_args
    argname=arglist{ii};
    found=0;
    % Is this name a flag? If so, set it
    if isfield(flagsreverse,argname)
      % Unset all other flags in this group
      flaggroup=defflags.(flagsreverse.(argname));
      for jj=1:numel(flaggroup)
        flags.(['do_',flaggroup{jj}])=0;
      end;
      
      flags.(flagsreverse.(argname))=argname;
      flags.(['do_',argname])=1;
      found=1;
    end;
      
    % Is this name the key of a key/value pair? If so, set the value.
    if isfield(defkeyvals,argname)      
      keyvals.(argname)=arglist{ii+1};
      ii=ii+1;
      found=1;
    end;
    
    if found==0
      error('%s: Unknown parameter: %s',callfun,argname);
    end;

    ii=ii+1;
  end;
