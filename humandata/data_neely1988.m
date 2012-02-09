function [F, tau]  = data_neely1988(varargin)
%DATA_NEELY1988 ABR wave V latency as function of center frequency
%   Usage: data = data_neely2010(flag)
%
%   Output parameters:
%     tau  :  Wave V latency
%     F    :  Center frequency of stimulus (same as input) 
%
%   `data_neely1988(flag)` returns data points based on equation 1 from Neely et
%   al. (1988). The figure produced is Neely et al. (1988)'s data fitted model and can be
%   compared to the raw data in figure 1.
%
%   The flag may be one of:
%
%     'noplot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%
%     'L'       Choose stimulus levels given as an array. Default value is `40:10:100`.
%
%     'F'       Choose stimulus center frequencies given as an
%               array. Defalt value is `1000:8000`.
%
%   Examples:
%   ---------
%
%   Figure can be displayed using:::
%
%     data_neely1988('plot');
%
%   References: neely1988latency

% Define input flags
definput.keyvals.L=40:10:100;
definput.keyvals.F=1000:8000; % Stimulus center frequency, Hz
definput.flags.plot = {'noplot','plot'};

[flags,kv]  = ltfatarghelper({},definput,varargin);
L       = kv.L/100; % Stimulus level, dB SPL

for i = 1:length(kv.L)
    tau(i,:)     = 5 + 12.9 * 5.^(-L(i)) * (kv.F/1000).^(-.413);
end
if flags.do_plot
    figure, semilogx(kv.F,tau','k-'),xlabel('CF'), ylabel('Latency [ms]')
end