function op1=amt_version(varargin)
%amt_version Help on the AMToolbox
%   Usage:  amt_version;
%           v=amt_version('version');
%           mlist=amt_version('modules');
%
%   `amt_version` displays some general help on the AMT.
%
%   `amt_version('version')` returns the version number.
%
%   `amt_version('modules')` returns a cell array of installed modules and
%   corresponding version numbers.
%
%   See also:  amt_start

%   AUTHOR : Peter Søndergaard. 
%   MODIFICATIONS: 09.07.2017 Piotr Majdak

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
  modinfo=amt_version('modules');
  for ii=1:length(modinfo);
    s=sprintf(' %-15s %7s  %s',modinfo{ii}.name,modinfo{ii}.version, ...
	      modinfo{ii}.description);
    amt_disp(s);
  end;

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
