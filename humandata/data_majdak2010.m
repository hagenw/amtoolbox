function data = data_majdak2010(varargin)
%DATA_MAJDAK2010 Listener specific localization performance
%   Usage: data = data_majdak2010(condition)
%
%   Output parameters:
%     data.id    : listener ID
%     data.mtx   : experimental data matrix conaining 9 colums
%                  col 1: target azimuth
%                  col 2: target elevation
%                  col 3: response azimuth
%                  col 4: response elevation
%                  col 5: lateral angle of target
%                  col 6: polar angle of target
%                  col 7: lateral angle of response
%                  col 8: polar angle of response
%
%   `data_majdak2010(condition)` returns listener-specific experimental data
%   from Majdak et al. (2010) testing localization performance for various
%   experimental methods.
% 
%   The *condition* flag may be one of:
%
%     'HMD_M'   Head-mounted display and manual pointing. Testing of naive
%               subjects.
%     'HMD_H'   Head-mounted display, head pointing, naive subjects.
%     'Dark_M'  Dark room, manual pointing, naive subjects.
%     'Dark_H'  Dark room, head pointing, naive subjects.
%     'Learn_M' Acoustic learning condition with manual pointing. This is the default.
%     'Learn_H' Acoustic learning condition with head pointing.
%
%   References: majdak2010methods

% AUTHOR: Robert Baumgartner

%% Check input options

% Define input flags
definput.flags.condition = {'Learn_M','Learn_H','HMD_M','HMD_H','Dark_M','Dark_H'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Extract data
if not(exist([mfilename '.mat'],'file'))
  disp(['Downloading ' mfilename ' from http://www.kfs.oeaw.ac.at/']);
  targetfn = fullfile(amtbasepath,'humandata',[mfilename '.mat']);
  sourcefn = ['http://www.kfs.oeaw.ac.at/research/experimental_audiology/projects/amt/' mfilename '.mat'];
  urlwrite(sourcefn,targetfn);
end
load(mfilename)

C = find(ismember(condition,flags.condition));

for ll = 1:length(subject)
  
  if not(isempty(subject(ll).expData{C}))
    data(ll).mtx = subject(ll).expData{C}(:,1:8);
  end
  data(ll).id = subject(ll).id;

end