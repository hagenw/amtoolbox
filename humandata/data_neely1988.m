function [F, tau]  = data_neely1988(varargin)
%DATA_ELBERLING2010 ABR wave V data as functon of level and sweeping rate
%   Usage: data = data_elberling2010(flag)
%
%   Output parameters:
%   tau  :  Wave V latency
%   F    :  Center frequency of stimulus (same as input) 
%
%   `data_neely1988(flag)´ returns data points based on equation 1 from Neely et
%   al. (1988). The figure produced is Neely et al. (1988)'s data fitted model and can be
%   compared to the raw data in figure 1 .
%
%   The flag may be one of:
%
%     'noplot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%
%     'L'       Choose stimulus levels - given as array
%
%     'F'       Choose stimulus center frequencies - given as array
%
%   Examples:
%   ---------
%
%   Figure can be displayed using:::
%
%     data_neely1988('plot','F', [1000 2000 5000 8000],'L', [40 60 80 100]);
%
%   References: Neeley et al (1988)

% Define input flags
definput.keyvals.L=40:10:100;
definput.keyvals.F=1000:8000; % Stimulus center frequency, Hz
definput.flags.plot = {'noplot','plot'};

[flags,kv, F, L]  = ltfatarghelper({'F', 'L'},definput,varargin);
L       = kv.L/100; % Stimulus level, dB SPL

for i = 1:length(kv.L)
    tau(i,:)     = 5 + 12.9 * 5.^(-L(i)) * (kv.F/1000).^(-.413);
end
if flags.do_plot
    figure, semilogx(kv.F,tau','k-'),xlabel('CF'), ylabel('Latency [ms]')
end