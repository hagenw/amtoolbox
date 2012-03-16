function tau  = data_neely1988(F,L,varargin)
%DATA_NEELY1988 ABR wave V data as functon of level and sweeping rate
%   Usage: tau = data_neely1988(flag)
%
%   Input parameters:
%     F    : Centre frequencies of stimulus
%     L    : Levels at centre frequencies
%
%   Output parameters:
%     tau  :  Wave V latency
%
%   `data_neely1988(F,L)` returns data points based on equation 1 from
%   Neely et al. (1988) where *F* is the centre frequencies and *L* are
%   the associated levels. 
%
%   The flag may be one of:
%
%     'noplot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%
%   Examples:
%   ---------
%
%   Figure XXX in Neely et al. (1988) can be reproduced using:::
%
%     F=[1000 2000 5000 8000];
%     tau=data_neely1988(F,[40 60 80 100]);
%     semilogx(F,tau','k-');
%     xlabel('CF');
%     ylabel('Latency [ms]')
%
%   References: neely1988latency

% Define input flags
%definput.keyvals.L=40:10:100;
%definput.keyvals.F=1000:8000; % Stimulus center frequency, Hz
%definput.flags.plot = {'noplot','plot'};

%[flags,kv]  = ltfatarghelper({},definput,varargin);

L = L/100; % Stimulus level, dB SPL

for ii = 1:length(L)
  tau(ii,:) = 5 + 12.9 * 5.^(-L(ii)) * (F/1000).^(-.413);
end
