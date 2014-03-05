function amtstart()
%AMTSTART   Start the Auditory Modelling Toolbox
%   Usage:  amtstart;
%
%   `amtstart` starts the Auditory Modelling toolbox. This command must be
%   run before using any of the function in the toolbox.
%
%   The auditory modelling toolbox depends on the Linear Time Frequency
%   Analysis Toolbox (LTFAT) to properly function. Therefore, you must issue
%   the `ltfatstart` command before you start AMT.
%
%   To configure default options for functions, you can use the
%   `ltfatsetdefaults` function in your startup script. A typical startup
%   file could look like::
%
%     addpath('/path/to/my/work/ltfat');
%     addpath('/path/to/my/work/amtoolbox');
%     ltfatstart;
%     amtstart;
%     ltfatsetdefaults('audspecgram','classic');
%
%   The last command wil configure |audspecgram| to display a classic
%   auditory spectrogram by default.
%
%   See also:  amthelp
%
%   References: ltfatnote015
  
%   AUTHOR : Peter L. Søndergaard.  

%% LTFAT package

% Search for LTAFT package
basepath=which('amtstart');
basepath=basepath(1:end-11);
if ~exist('ltfatstart','file')
  ltfatpath=fullfile(basepath,'thirdparty','LTFAT');
  if exist(ltfatpath,'dir')
    addpath(ltfatpath);
  end
end

% Start LTFAT
disp('*** Starting LTFAT ***');
if exist('ltfatstart','file')
  ltfatstart;
else
  error(['LTFAT package could not be found. Unable to continue.' 10 ...
        'Download LTFAT from http://ltfat.sourceforge.net ' 10 ...
        'and copy to amtoolbox/thirdparty/LTFAT.']); 
end

% Check for the correct version. 
s=ltfathelp('version'); 
s_r='1.0.9'; % set the required version
v=sscanf(s,'%d.%d.%d'); v(4)=0;
v_r=sscanf(s_r,'%d.%d.%d');
if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
    error(['You need LTFAT >= ' s_r ' to work with AMT. ' ...
      'Please update your package from http://ltfat.sourceforge.net ']);
end  
     
%% SOFA package

% Search for SOFA package
basepath=which('amtstart');
basepath=basepath(1:end-11);
if ~exist('SOFAstart','file')
  sofapath=fullfile(basepath,'thirdparty','SOFA','API_MO');
  if exist(sofapath,'dir')
    addpath(sofapath);
  end
end

% Start SOFA
disp('*** Starting SOFA ***');
if exist('SOFAstart','file')
  SOFAdbPath(fullfile(basepath,'hrtf'));
  SOFAdbURL('http://www.sofacoustics.org/data/amt');
  SOFAstart;
else
  disp(['SOFA package could not be found. Continue without SOFA support.']);
  disp(['For SOFA support please download the package ' ...
        'from http://sofacoustics.sourceforge.net ' ...
        'and copy to amtoolbox/thirdparty/SOFA.']); 
end

%% SFS package

% Search for the package
basepath=which('amtstart');
basepath=basepath(1:end-11);
if ~exist('SFS_start','file')
  sfspath=fullfile(basepath,'thirdparty','SFS');
  if exist(sfspath,'dir')
    addpath(sfspath);
  end
end

% Start 
disp('*** Starting SFS ***');
if exist('SFS_start','file')
  SFS_start;
  s=SFS_version; s_r='0.2.4'; % set the required version
  disp(['SFS, version ' s]);
  v=sscanf(s,'%d.%d.%d'); v(4)=0;
  v_r=sscanf(s_r,'%d.%d.%d');
  if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
      error(['You need SFS >= ' s_r ' to work with AMT. ' ...
        'Please update your package from https://github.com/sfstoolbox/sfs ']);
  end  
else
  disp(['SFS package could not be found. Continue without SFS support.']);
  disp(['For SFS support please download the package ' ...
        'from https://github.com/sfstoolbox/sfs ' ...
        'and copy to amtoolbox/thirdparty/SFS.']); 
end

%% Start AMT
disp('*** Starting AMT ***');  
% --- general settings ---
% Print the banner at startup?
printbanner=1;

% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('amtstart');
basepath=basepath(1:end-11);
if exist('addpath','var')>0
  addpath(basepath);
else
  path(path,basepath);
end
bp=[basepath,filesep];

% Load the version number
[FID, MSG] = fopen ([bp,'amtoolbox_version'],'r');
if FID == -1
    error(MSG);
else
    amt_version = fgetl (FID);
    fclose(FID);
end

% -----------  install the modules -----------------

modules={};
nplug=0;

% List all files in base directory
d=dir(basepath);

for ii=1:length(d)
  if d(ii).isdir
    if ~(d(ii).name(1)=='.')

      name=d(ii).name;
      
      % The file is a directory and it does not start with '.' This could
      % be a module      
      if exist([bp,name,filesep,name,'init.m'],'file')
	% Set 'status' to zero if the module forgets to define it.
	status=0;
	module_version=amt_version;
        addpath([bp,name]);

	eval([name,'init']);
        if status>0
          if status==1
            nplug=nplug+1;
            modules{nplug}.name=name;
            modules{nplug}.version=module_version;
          end;
	else
	  rmpath([bp,name]);
	end;
      end;	

    end;
  end;
end;

% Check if Octave was called using 'silent'
%if isoctave
%  args=argv;
%  for ii=1:numel(args)
%    s=lower(args{ii});
%    if strcmp(s,'--silent') || strcmp(s,'-q')
%      printbanner=0;
%    end;
%  end;
%end;

if printbanner
  disp(['AMT version ',amt_version,'. Copyright 2009 - 2014 Peter L. Søndergaard. For help, please type "amthelp".'])
end;



%% ---------- load information into ltfathelp ------------

% As comp is now in the path, we can call ltfatarghelper
ltfatsetdefaults('amthelp','versiondata',amt_version,...
                 'modulesdata',modules);
