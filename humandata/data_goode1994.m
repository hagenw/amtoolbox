function data = data_goode1994(varargin)
%DATA_GOODE1994 Returns data points from the Goode et al. (1994) paper
%   Usage: data = data_goode1994_data(flag)
%
%   DATA_GOODE1994(flag) returns data points from the paper by Goode et
%   al. (1994). Currently, only Figure 1 with condition 104 dB is
%   supported.
%
%   The flag may be one of:
%
%-    'noplot'       - don't plot, only return data. This is the default.
%
%-    'plot'         - plot the data.
%
%-    'fig1_104'     - return data from Fig. 1. for 104 dB SPL, stapes
%                      footplate diplacement at 104 dB SPL. This is
%                      the default (and currently only option).
%
%R  goode1994nkf

%   AUTHOR: Peter L. Soendergaard


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type={'fig1_104'};
definput.flags.plot = {'noplot','plot'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%% ------ Data points from the paper ------------------------------------
%
% Data for the given figure

if flags.do_fig1_104
  % data derived from Goode et al. 1994, Figure 1
  data = [...
    400	0.19953; ...
    600	0.22909; ...
    800	0.21878; ...
    1000 0.15136; ...
    1200 0.10000; ...
    1400 0.07943; ...
    1600 0.05754; ...
    1800 0.04365; ...
    2000 0.03311; ...
    2200 0.02754; ...
    2400 0.02188; ...
    2600 0.01820; ...
    2800 0.01445; ...
    3000 0.01259; ...
    3500 0.00900; ...
    4000 0.00700; ...
    4500 0.00457; ...
    5000 0.00500; ...
    5500 0.00400; ...
    6000 0.00300; ...
    6500 0.00275];

  if flags.do_plot
    figure;
    loglog(data(:,1)/1000,data(:,2),'+');
    xlabel('Frequency (kHz)');
    ylabel('Peak-to-peak displacement ({\mu}m)');
    %axis([0.1,10,min(data(:,2))*1e6,max(data(:,2))*1e6]);
    axis([0.1,10,0.001,10]);
  end;
end;

