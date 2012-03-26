function plotroenne2012_tonebursts(waveVlat,click_latency)
%PLOTROENNE2012_TONEBURSTS plots Rønne et al. (2012) Fig. 5
%   Usage: plotroenne2012_tonebursts(flag)
%
%   `plotroenne2012_tonebursts(waveVlat,click_latency)` plots the output
%   from |roenne2012_tonebursts|_ in a similar way as Fig. 5 from the Rønne
%   et al. (2012) ABR model.
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%  
%   References: roenne2012modeling neely1988latency harte2009comparison zilany2007representation

if 1
  waveVlat = [10.7666666666667,9.96666666666667, ...
              9.80000000000000,8.96666666666667,8.03333333333333,8, ...
              7.93333333333333;...
              
              10.0333333333333,9.46666666666667, ...
              8.86666666666667,8.60000000000000,8.50000000000000, ...
              7.93333333333333,7.86666666666667;...
              9.60000000000000, ...
                      9.06666666666667,8.66666666666667,8.30000000000000, ...
                      7.86666666666667,7.73333333333333,7.36666666666667; ...
                      8.70000000000000,8.26666666666667,7.90000000000000, ...
                      7.60000000000000,7.33333333333333,7.10000000000000, ...
                      6.86666666666667;7.23333333333333,6.96666666666667, ...
                      6.76666666666667,6.56666666666667,6.40000000000000, ...
                      6.26666666666667,6.16666666666667;6.83333333333333, ...
                      6.60000000000000,6.43333333333333,6.30000000000000, ...
                      6.20000000000000,6.10000000000000,6.03333333333333];
  
  click_latency    = [6.76666666666667,6.53333333333334, ...
                      6.30000000000000,6.13333333333333,6,5.90000000000000, ...
                      5.86666666666666];
  
  %[flags,kv,waveVlat,click_latency]  = ...
  %    ltfatarghelper({'waveVlat','click_latency'},definput,varargin);
  
end;
  
% sampling frequency of simulations
fs = 30e3;

% Center Frequencies of stimuli
CF = [1000, 1500, 2000, 3000 ,6000, 8000];                 

%% Plot

% plot click latencies
hold all;
semilogx(10e3, click_latency(1),'ko','linewidth',2);
semilogx(10e3, click_latency(2),'kv','linewidth',2);
semilogx(10e3, click_latency(3),'kd','linewidth',2);
semilogx(10e3, click_latency(4),'k*','linewidth',2);
semilogx(10e3, click_latency(5),'kx','linewidth',2);
semilogx(10e3, click_latency(6),'kp','linewidth',2);
semilogx(10e3, click_latency(7),'k^','linewidth',2);
legend('40 dB','50 dB','60 dB','70 dB','80 dB','90 dB','100 dB');

% plot toneburst latencies
set(gca,'fontsize',12);
hold all;
semilogx(CF, waveVlat(:,1),'-ko','linewidth',2);
semilogx(CF, waveVlat(:,2),'-kv','linewidth',2);
semilogx(CF, waveVlat(:,3),'-kd','linewidth',2);
semilogx(CF, waveVlat(:,4),'-k*','linewidth',2);
semilogx(CF, waveVlat(:,5),'-kx','linewidth',2);
semilogx(CF, waveVlat(:,6),'-kp','linewidth',2);
semilogx(CF, waveVlat(:,7),'-k^','linewidth',2);

% Plot Neely et al (1988) data
F=1e3:8e3;
tau = data_neely1988(F,40:10:100);
hold on;
semilogx(F,tau','k--');
text(10e3,8, 'Click','fontsize',12);

hold on; 

plot(12e3, [5.9, 6.6, 7.6],'ok','MarkerFaceColor','k');
text(13e3,5.88, '95.2');
text(13e3,6.59, '75.2');
text(13e3,7.6, '55.2');

% Plot Harte et al (2009) data
F       = 800:10e3;
tauHarte =5+ 11.09*(F/1000).^(-0.37)/2;
semilogx(F,tauHarte,'k:','linewidth',1.5);

% Plot init
box on;
ylabel('Latency of wave V [ms]');
xlabel('Frequency of toneburst [kHz]');
set(gca,'fontsize',12);
set(gca,'XTick',CF,'XTickLabel',CF/1000);
set(gca,'YTick',1:15,'YTickLabel',1:15);
axis([800 16000 5.5 12.5]);

