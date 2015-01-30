function data = data_majdak2013(varargin)
%DATA_MAJDAK2013 Listener specific localization in saggital planes
%   Usage: data = data_majdak2013(condition)
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
%   `data_majdak2013(condition)` returns listener-specific experimental data
%   from Majdak et al.  (2013) testing localization performance in sagittal
%   planes for low-pass filtered and spectrally warped DTFs.
%
%   The *condition* flag may be one of:
%
%     'BB'   Broadband DTFs (baseline condition). This is the default.
%     'LP'   Low-pass filtered (at 8.5kHz) DTFs
%     'W'    Spectrally warped (2.8-16kHz warped to 2.8-8.5kHz) DTFs
%
%   References: majdak2013spatstrat

% AUTHOR: Robert Baumgartner

%% Check input options

% Define input flags
definput.flags.condition = {'BB','LP','W'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Extract data
% if not(exist([mfilename '.mat'],'file'))
%   amtdisp(['Downloading ' mfilename ' from http://www.kfs.oeaw.ac.at/']);
%   targetfn = fullfile(amtbasepath,'humandata',[mfilename '.mat']);
%   sourcefn = ['http://www.kfs.oeaw.ac.at/research/experimental_audiology/projects/amt/' mfilename '.mat'];
%   urlwrite(sourcefn,targetfn);
% end
% load(mfilename)
x=amtload('majdak2013','data.mat');

C = find(ismember(x.condition,flags.condition));

for ll = 1:length(x.subject)
  
  data(ll).mtx = x.subject(ll).expData{C}(:,1:8);
  data(ll).id = x.subject(ll).id;

end