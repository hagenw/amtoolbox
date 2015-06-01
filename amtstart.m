function amtstart(varargin)
%AMTSTART   Start the Auditory Modeling Toolbox
%   Usage:  amtstart;
%
%   `amtstart` starts the Auditory Modeling Toolbox. This command must be
%   run before using any of the function in the toolbox.
%
%   The AMT depends on the Linear Time Frequency Analysis Toolbox (LTFAT). 
%   You must first download LTFAT from
%   http://ltfat.sourceforge.net/ and unpack the downloaded file. 
%   In the AMT, there is a pre-prepared directory `thirdparty/ltfat`
%   where the LTFAT can be stored. Alternatively, set the path to your
%   LTFAT installation to the search path of Matlab/Octave.
%
%   In order to run all the AMT functionality, you will need to:
%   
%   1) install SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. `thirdparty/SOFA`)
%   2) install SFS Toolbox from https://github.com/sfstoolbox/sfs
%   4) install Python >2.6 with `numpy` and `scipi` packages. On Linux, use `sudo apt-get install python-scipy python-numpy`
%   3) run `amtmex` and compile successfully
%   5) run `make` (Linux) or `make.bat` (Windows) in `src/verhulst`
%   6) have the Optimization Toolbox for Matlab installed
% 
%   Some of the AMT functions require a large processing time. Depending on the machine and the model, it might take even days. Thus, some AMT functions provide caching of calculated results. If you don't want to wait and just take a look at the results: download the cached data from https://sourceforge.net/projects/amtoolbox/files/, unzip to the root AMT directory, and run the particular AMT function.
%
%   Most of the models require auxiliary data. The AMT will download these data on-demand. 
%   Some of the models require HRTFs. The AMT will download alse these HRTFs on-demand.
%   If you want to run the AMT completely offline, download the auxiliary data and HRTFs ahead. 
%   The download URL for the auxiliary data is given by `amtauxdataurl`. 
%   The target directory for the auxiliary data is given by `amtauxdatapath`. 
%   The download URL for the HRTFs is given by `SOFAdbURL`.
%   The target directory for the HRTFs is given by `SOFAdbPath`. 
%
%   `amtstart('documentation')` starts the AMT in the documentation compiling
%   mode. The progress output will be suppressed.
%
%   `amtstart('silent')` starts the AMT in the silent mode where all output but figures will be suppressed.
%
%   `amtstart('verbose')` starts the AMT in the verbose mode and all output will be dispplayed. This is the default mode. 
% 
%   See also:  amtmex amtflags amtload amtcache
%
  
%   AUTHOR : Peter L. Soendergaard, Piotr Majdak 


%% Start AMT
bp=amtbasepath;

% Load the version number
[FID, MSG] = fopen ([bp,'amtoolbox_version'],'r');
if FID == -1
    error(MSG);
else
    amt_version = fgetl (FID);
    fclose(FID);
end

% Check if 'silent' present in the flags
silent=0;
if isoctave, args=argv; else args=varargin; end
 for ii=1:numel(args)
   s=lower(args{ii});
   if strcmp(s,'silent') || strcmp(s,'-q')
     silent=1;
   end;
 end;
% end;


if ~silent
  disp('  ');
  disp(['AMT version ',amt_version,'. (C) Peter L. Soendergaard and Piotr Majdak.']);
  disp('  ');
  disp('Starting toolboxes...');
end;

%% LTFAT package

% Search for LTAFT package
basepath=which('amtstart');
basepath=basepath(1:end-11);
if ~exist('ltfatstart','file')
  ltfatpath=fullfile(basepath,'thirdparty','ltfat');
  if exist(ltfatpath,'dir')
    addpath(ltfatpath);
  end
end

% Start LTFAT
% if ~silent, disp('*** Starting LTFAT ***'); end
if exist('ltfatstart','file')
  if silent, ltfatstart(0); else ltfatstart; end;
else
  error(['LTFAT package could not be found. Unable to continue.' 10 ...
        'Download LTFAT from http://ltfat.sourceforge.net ' 10 ...
        'and copy to amtoolbox/thirdparty/ltfat.']); 
end

% Check for the correct version. 
s=ltfathelp('version'); 
s_r='2.0.0'; % set the required version
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
% if ~silent, disp('*** Starting SOFA ***'); end
if exist('SOFAstart','file')
  SOFAdbPath(fullfile(basepath,'hrtf'));
  SOFAdbURL('http://www.sofacoustics.org/data/amt');
  if silent, SOFAstart('silent'); else SOFAstart; end
	warning('off','SOFA:upgrade');	% disable warning on upgrading older SOFA files
	warning('off','SOFA:load'); % disable warnings on loading SOFA files
else
  if ~silent,
  disp('SOFA package could not be found. Continue without SOFA support.');
  disp(['For SOFA support please download the package ' ...
        'from http://sofacoustics.sourceforge.net ' ...
        'and copy to amtoolbox/thirdparty/SOFA.']); 
  end
end

%% SFS package

% Search for the package
basepath=which('amtstart');
basepath=basepath(1:end-11);
if ~exist('SFS_start','file')
  sfspath=fullfile(basepath,'thirdparty','sfs');
  if exist(sfspath,'dir')
    addpath(sfspath);
  end
end

% Delete rms.m from the SFS package because of naming conflict
sfspath=fileparts(which('SFS_start.m'));
if exist(fullfile(sfspath,'SFS_general','rms.m'),'file'),
	delete(fullfile(sfspath,'SFS_general','rms.m'));
end

% Start 
% if ~silent, disp('*** Starting SFS ***'); end
if exist('SFS_start','file')
  SFS_start;
  s=SFS_version; s_r='1.0.0'; % set the required version
  if ~silent, disp(['Sound Field Synthesis Toolbox, version ' s]); end
  v=sscanf(s,'%d.%d.%d'); v(4)=0;
  v_r=sscanf(s_r,'%d.%d.%d');
  if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
      error(['You need SFS >= ' s_r ' to work with AMT. ' ...
        'Please update your package from https://github.com/sfstoolbox/sfs ']);
  end  
	
elseif ~silent, 
  disp('SFS package could not be found. Continue without SFS support.');
  disp(['For SFS support please download the package ' ...
        'from https://github.com/sfstoolbox/sfs ' ...
        'and copy to amtoolbox/thirdparty/sfs.']); 
end

%% Install AMT modules
% A directory called DIRNAME containing a file 'DIRNAMEinit.m' is
% considered as a module. 
% DIRNAMEinit.m must set the variable 'status' with the following value:
%  0: disabled module, don't add to the search path
%  >0: add to the search path.

% add root of the AMT to the path
addpath(basepath);

modules={};
nplug=0;

% List all files in base directory
d=dir(basepath);

for ii=1:length(d)
  if d(ii).isdir
    if ~(d(ii).name(1)=='.')      
      % The file is a directory and it does not start with '.' This could be a module      
      name=d(ii).name;
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

%% define default start-up behaviour
flags=amtflags(varargin); % amtdisp and other amt-related functions work now!

%% ---------- load information into ltfathelp ------------

% As comp is now in the path, we can call ltfatarghelper
ltfatsetdefaults('amthelp','versiondata',amt_version,...
                 'modulesdata',modules);

%% Initialize aux data, cache, and display starting information
amtdisp('  ');
% amtdisp('  ');
amtdisp('*** AMT ready to go! ***'); 
amtdisp(['Auxiliary data (local): ' amtauxdatapath]);
amtdisp(['Auxiliary data (web): ' amtauxdataurl]);
if strcmp(flags.cachemode,'global'), flags.cachemode='normal'; end
amtcache('setMode',flags.cachemode);
switch flags.cachemode
  case 'normal'
    amtdisp('Cache mode: Download precalculated results.');
    amtdisp('            exp_model(...)          shows precalculated results');
    amtdisp('            exp_model(...,''redo'') enforces recalculation');
  case 'localonly'
    amtdisp('Cache mode: Use local cache or recalculate. Do not connect to remote cache.');
  case 'cached'
    amtdisp('Cache mode: Enforce using cache. Do not recalcalculate.');
  case 'redo'
    amtdisp('Cache mode: Recalculate always (be patient!).');
end

