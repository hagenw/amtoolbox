function [acrossSubdSRTs,std_acrossSubSRTs]  = data_joergensen2011(varargin)
%DATA_JOERGENSEN2011 XXX
%   Usage: ur = data_joergensen2011()
%
%   `[acrossSubdSRTs,std_acrossSubSRTs]=data_joergensen2011('figXXX`)
%   returns the data from Fig. XXX ...
%
%   Examples:
%   ---------
%
%   XXX Please provide a plot here to visualize the relevant figure:::
%
%     [ur,fs]  = data_joergensen2011;
%     plot((0:length(ur)-1)/fs,ur);
%     xlabel('Time / seconds');
%     ylabel('Amplitude');
%
%   References: joergensen2011predicting
  

acrossSubdSRTs = ...
    [0, 1.3750, 1.3750, 1.6250, 1.8750, 2.7083];

std_acrossSubSRTs = ...

    [0.4513, 0.7500, 0.2500, 0.4590, 0.3436, 0.5159];