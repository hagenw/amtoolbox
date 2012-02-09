function [delay,data_mean,data_std]  = data_elberling2010(varargin)
%DATA_ELBERLING2010 ABR wave V data as functon of level and sweeping rate
%   Usage: data = data_elberling2010(flag)
%
%   Output parameters:
%     delay      : "x-axis" - sweeping rate delay between 710Hz and 5700Hz
%     data_mean  : Mean of data 
%     data_std   : Standard deviation of data
%
%   `data_elberling2010(flag)` returns data points from Elberling et
%   al. (2010)
%
%   The flag may be one of:
%
%     'noplot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%  
%     'fig4'    Data from Fig. 4, Amplitude of wave V. This is the default.
%
%     'fig5'    Data from Fig. 5, Latency of wave V
%
%   Examples:
%   ---------
%
%   Figure 4 can be displayed using:::
%
%     data_elberling2010('fig4','plot');
%
%   Figure 5 can be displayed using:::
%
%     data_elberling2010('fig5','plot');
%
%   References: elberling2010evaluating
    
% Define input flags
definput.flags.type={'fig4','fig5'};
definput.flags.plot = {'noplot','plot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

ElberlingDelay = [0 1.86 2.56 3.32 4.14 5.04]; % Delays form Elberling et al (2010) "x-axis"

ElberlingDelay2 = [ElberlingDelay;ElberlingDelay;ElberlingDelay];

% Set the fontsize and color.
ftz = 12;
%color = [0.7,0.7,0.7];
color = [0,0,0];
fs=3e4;
if flags.do_fig4
  
  ElberlingMean    = [256 311 339 349 382 391; ...
                      368 560 613 645 631 581 ; ...
                      407 624 654 596 516 389]';
  
  ElberlingStdDev  = [61   86  86  72  83  91; ...
                      92  135 119  79  79 100 ; ...
                      91  162 118 118 96  120]';
  
  if flags.do_plot
    figure;

    errorbar(ElberlingDelay2',ElberlingMean,ElberlingStdDev/sqrt(20),'-',...
             'linewidth',1.5,'color',color)

    %set(gca,'fontsize',12);
    
    %ylim([0 800]);
    axis([-1.2 5.5 0 800]);
    xlabel('Change of delay [ms]');
    ylabel('ABR amplitude [nv]');
    text(-.7,ElberlingMean(1,1), '20','fontsize',ftz)
    text(-.7,ElberlingMean(1,2), '40','fontsize',ftz)
    text(-.7,ElberlingMean(1,3), '60','fontsize',ftz)
    text(-.9,ElberlingMean(1,3)+70, 'dB nHL','fontsize',ftz)
    text(-.2,50, 'Click','fontsize',ftz)
    text(ElberlingDelay(2),75, '1','fontsize',ftz)
    text(ElberlingDelay(3),75, '2','fontsize',ftz)
    text(ElberlingDelay(4),75, '3','fontsize',ftz)
    text(ElberlingDelay(5),75, '4','fontsize',ftz)
    text(ElberlingDelay(6),75, '5','fontsize',ftz)
    text(3,40, 'Chirps','fontsize',ftz)
    

  end;
  
  data_mean = ElberlingMean;
  data_std  = ElberlingStdDev;
end;

if flags.do_fig5
  
  ElberlingLatencyMean = [7.60 7.23 7.12 6.97 6.86 6.66; ...
                      6.59 6.15 5.96 5.76 5.48 4.98; ...
                      5.88 4.97 4.57 4.12 3.50 2.87]';
  
  ElberlingLatencyStdDev  = [0.79 0.36 0.32 0.34 0.36 0.36; ...
                      0.74 0.29 0.23 0.27 0.34 0.59; ...
                      0.25 0.28  0.35 0.51 0.51 0.38]';
  
  if flags.do_plot
    figure; 

    errorbar(ElberlingDelay2',ElberlingLatencyMean,ElberlingLatencyStdDev,'-',...
             'linewidth',1.5,'color',color)

    %set(gca,'fontsize',ftz);
    axis([-1.2 6.5 0 10]);
    xlabel('Change of delay [ms]');
    ylabel('ABR latency [ms]')
    text(-.7,5.88, '60','fontsize',ftz,'color',color)
    text(-.7,6.59, '40','fontsize',ftz,'color',color)
    text(-.7,7.6, '20','fontsize',ftz,'color', color)
    text(-.9,8.7, 'dB nHL','fontsize',ftz,'color',color)
    text(5.5,ElberlingLatencyMean(6,1)/fs*1000-15, '20','fontsize',ftz)
    text(5.5,ElberlingLatencyMean(6,2)/fs*1000-15, '40','fontsize',ftz)
    text(5.5,ElberlingLatencyMean(6,3)/fs*1000-15, '60','fontsize',ftz)
    text(5.3,ElberlingLatencyMean(6,1)/fs*1000-15+.6, 'dB nHL','fontsize',ftz)
    text(-.2,.5, 'Click','fontsize',ftz)
    text(ElberlingDelay(2),1, '1','fontsize',ftz)
    text(ElberlingDelay(3),1, '2','fontsize',ftz)
    text(ElberlingDelay(4),1, '3','fontsize',ftz)
    text(ElberlingDelay(5),1, '4','fontsize',ftz)
    text(ElberlingDelay(6),1, '5','fontsize',ftz)
    text(3,.5, 'Chirps','fontsize',ftz)
    %box on;
    set(gca,'fontsize',ftz);
    %ylim([0 800]);
  end;
  
  data_mean = ElberlingLatencyMean;
  data_std  = ElberlingLatencyStdDev;
end;

delay = ElberlingDelay;