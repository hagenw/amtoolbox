function [outx,outy]  = data_joergensen2011(varargin)
%DATA_JOERGENSEN2011 XXX
%   Usage:  = data_joergensen2011()
%
%   `[outx,outy]=data_joergensen2011('figXXX`)
%   returns data points from the Joergensen2011 paper. The flag may be one
%   of:
%
%     'noplot'        Don't plot, only return data. This is the default.
%
%     'plot'          Plot the data.
%
%     'figXXX'        Describe output data.
%
%     'figXXX'        Describe output data.
%
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
  
definput.flags.plot = {'noplot','plot'};
definput.flags.type = {'figXXX','figYYY'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_figXXX
    acrossSubdSRTs = ...
        [0, 1.3750, 1.3750, 1.6250, 1.8750, 2.7083];    
    std_acrossSubSRTs = ...        
        [0.4513, 0.7500, 0.2500, 0.4590, 0.3436, 0.5159];
    
    outx=acrossSubdSRTs;
    outy=std_acrossSubSRTs;
    
    if flags.do_plot
        % Put the visualization of the data here
        
    end;
    
end;

if flags.do_figYYY
    
    % repeat
    
end;

