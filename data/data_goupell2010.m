function data = data_goupell2010(varargin)
%DATA_GOUPELL2010 Localization performance in sagittal planes
%   Usage: data = data_goupell2010(condition)
%          data = data_goupell2010(lat, dlat, condition)
%
%   Listener-specific experimental data from Goupell et al. (2010) testing
%   localization performance in sagittal planes for various numbers of
%   channels of a GET vocoder.
%
%   The *condition* flag may be one of:
%
%     'BB'   Broadband DTFs (baseline condition). This is the default.
%     'CL'   Click trains with unlimited number of channels
%     'N24'  24 vocoder channels
%     'N18'  18 vocoder channels
%     'N12'  12 vocoder channels
%     'N9'   9 vocoder channels
%     'N6'   6 vocoder channels
%     'N3'   3 vocoder channels
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
%   References goupell2010numchan

% AUTHOR: Robert Baumgartner

%% Check input options

% Define input flags
definput.flags.condition = {'BB','CL','N24','N18','N12','N9','N6','N3'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);

%% Extract data
x=amt_load('goupell2010','data.mat');

C = find(ismember(x.condition,flags.condition));

for ll = 1:length(x.subject)
  
  data(ll).mtx = x.subject(ll).expData{C}(:,1:8);
  data(ll).id = x.subject(ll).id;

end