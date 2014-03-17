function data = data_majdak2013ctc(varargin)
%DATA_MAJDAK2013CTC Listener-specific localization in sagittal planes
%   Usage: data = data_majdak2013ctc(condition)
%          data = data_majdak2013ctc(lat, dlat, condition)
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
%   `data_majdak2013ctc(condition)` returns listener-specific experimental
%   data from Majdak et al. (2013) testing localization performance in
%   sagittal planes for repeated HRTF measurements motivated by CTC binaural
%   synthesis.
%
%   The *condition* flag may be one of:
%
%     'Learn' Last 300 trials of acoustical training with visual feedback.
%     'A'     First HRTF measurement. This is the default.
%     'B'     Second HRTF measurement.
%
%
%   References: majdak2013ctc

% AUTHOR: Robert Baumgartner

%% Check input options

% Define input flags
definput.flags.condition = {'A','B','Learn'};

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