function data = data_majdak2013(varargin)
%DATA_majdak2013 Listener-specific experimental data from Majdak et al.
%(2013) testing localization performance in sagittal planes for low-pass
%filtered and spectrally warped DTFs.
%
%   Usage: data = data_majdak2013(condition)
%          data = data_majdak2013(lat, dlat, condition)
%
%   The *condition* flag may be one of:
%
%     'BB'   Broadband DTFs (baseline condition). This is the default.
%     'LP'   Low-pass filtered (at 8.5kHz) DTFs
%     'W'    Spectrally warped (2.8-16kHz warped to 2.8-8.5kHz) DTFs
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
%   References: majdak2013spatstrat

% AUTHOR: Robert Baumgartner

%% Check input options

% Define input flags
definput.flags.condition = {'BB','LP','W'};

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
  
  data(ll).mtx = subject(ll).expData{C}(:,1:8);
  data(ll).id = subject(ll).id;

end