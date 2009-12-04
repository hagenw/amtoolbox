function op1=amthelp(q,modulename)
%AMTHELP Help on the AMT model
%   Usage:  amthelp;
%           v=amthelp('version');
%           mlist=amthelp('modules');
%
%   AMTHELP displays some general help on the AMT toolbox and model.
%
%   AMTHELP('version') returns the version number.
%
%   AMTHELP('modules') returns a cell array of installed modules and
%   corresponding version numbers.
%
%   AMTHELP('authors') lists the name of the authors.
%
%   See also:  amtstart

%   AUTHOR : Peter Soendergaard.  
%   TESTING: NA
%   REFERENCE: NA

global AMT_CONF;

% Verify that AMT_CONF has been initialized
if numel(AMT_CONF)==0
  disp(' ');
  disp('--- AMT - Auditory Modelling Toolbox. ---');
  disp(' ')
  disp('To start the toolbox, call AMTSTART as the first command.');
  disp(' ');
  return;
end;

if nargin==0
  disp(' ');
  disp('--- AMT - Auditory Modelling Toolbox. ---');
  disp(' ')

  disp(['Version ',AMT_CONF.amt_version]);
  disp(' ');
  disp('Installed modules:');
  disp(' ');
  disp('Name:            Version:  Description');
  modinfo=amthelp('modules');
  for ii=1:length(modinfo);
    s=sprintf(' %-15s %7s  %s',modinfo{ii}.name,modinfo{ii}.version, ...
	      modinfo{ii}.description);
    disp(s);
  end;

  disp(' ')
  if isoctave
    disp('Type amthelp("modulename") where "modulename" is the name of one');
    disp('of the modules to see help on that module.');
    
  else
    disp('Type "help modulename" where "modulename" is the name of one')
    disp('of the modules to see help on that module.') 

  end; 
  disp(' ');
  disp('For other questions, please don''t hesitate to send an email to amt-help@lists.sourceforge.net.'); 
    
  return
end;
  
if ischar(q);

  bp=AMT_CONF.basepath;

  switch(lower(q))
    case {'v','version'}
      if nargin==1
	op1=AMT_CONF.amt_version;
      else
	for ii=1:length(AMT_CONF.modules);
      % strcmpi does not exist in older versions of Octave,
      % therefore use strcmp(lower(...
	  if strcmp(lower(modulename),AMT_CONF.modules{ii}.name)
	    found=1;
	    op1=AMT_CONF.modules{ii}.version;
	    return
	  end;
	end;

	% If we get here, it means that no module matched 'modulename'.
	error('Unknown module name.');
				   
      end;
    case {'m','modules'}
      op1={};
      for ii=1:length(AMT_CONF.modules);

	p=AMT_CONF.modules{ii};

	% Get the first line of the help file
	[FID, MSG] = fopen ([bp,p.name,filesep,'Contents.m'],'r');
        if FID==-1
          error('Module %s does not contain a Contents.m file.',p.name);
	end;
	firstline = fgetl (FID);
	fclose(FID);

    
	% Load the information into the cell array.	
	op1{ii}.name=p.name;
	op1{ii}.version=p.version;
	op1{ii}.description=firstline(2:end);
        
      end;
    case {'a','authors'}
      disp('Peter L. Soendergaard, Centre for Applied Hearing Research,');
      disp('                       Technical University of Denmark.');
      disp(' ');
      disp('Tobias May,            Phillips Research Laboratories,');
      disp('                       Eindhoven.');
      
    otherwise
      found=0;
      for ii=1:length(AMT_CONF.modules);
	if strcmp(lower(q),AMT_CONF.modules{ii}.name)
	  found=1;
	  p=AMT_CONF.modules{ii};

	  if isoctave
	    s=pwd;
	    cd([bp,p.name]);
	    help Contents
	    cd(s);
	  else
	    help(lower(q))
	  end;
	end;
      end;  
      if ~found
	error('Unknown command or module name.');
      end;
  end;
end;
