function op1=amt_help(varargin)
%amt_help Help on the AMToolbox
%   Usage:  amt_help;
%           v=amt_help('version');
%           mlist=amt_help('modules');
%
%   `amt_help` displays some general help on the AMToolbox.
%
%   `amt_help('version')` returns the version number.
%
%   `amt_help('modules')` returns a cell array of installed modules and
%   corresponding version numbers.
%
%   See also:  amt_start

%   AUTHOR : Peter Søndergaard.  
%   TESTING: NA
%   REFERENCE: NA

% Verify that LTFAT has been installed
if ~exist('ltfatarghelper','file')  
    amt_disp('');
    amt_disp('--- AMTOOLBOX - The Auditory Modeling toolbox. ---');
    amt_disp('')
    error(['The toolbox require the LTFAT toolbox to properly function. ' ...
         'Please download and install it from http://ltfat.sourceforge.net,' ...
         'and then call the LTFATSTART command BEFORE you use the AMT.'])
end;

bp=amt_basepath;

definput.keyvals.versiondata=[];
definput.keyvals.modulesdata=[];
definput.flags.mode={'general','version','modules','authors'};

[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_general
  amt_disp(' ');
  amt_disp('--- AMT - The Auditory Modeling Toolbox. ---');
  amt_disp(' ')

  amt_disp(['Version ',kv.versiondata]);
  amt_disp(' ');
  amt_disp('Installed modules:');
  amt_disp(' ');
  amt_disp('Name:            Version:  Description');
  modinfo=amt_help('modules');
  for ii=1:length(modinfo);
    s=sprintf(' %-15s %7s  %s',modinfo{ii}.name,modinfo{ii}.version, ...
	      modinfo{ii}.description);
    amt_disp(s);
  end;

  amt_disp(' ')
  if isoctave
    amt_disp('Type amt_help("modulename") where "modulename" is the name of one');
    amt_disp('of the modules to see help on that module.');
    
  else
    amt_disp('Type "help modulename" where "modulename" is the name of one')
    amt_disp('of the modules to see help on that module.') 

  end; 
  amt_disp(' ');
  amt_disp('For other questions, please don''t hesitate to send an email to amtoolbox-help@lists.sourceforge.net.'); 
    
end;
  
if flags.do_version
  op1=kv.versiondata;
end;

if flags.do_modules
  op1={};
  for ii=1:numel(kv.modulesdata)
    
    p=kv.modulesdata{ii};
    
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
end;
