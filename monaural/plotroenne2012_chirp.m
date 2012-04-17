function plotroenne2012_chirp(waveVamp, waveVlat, varargin)
%PLOTROENNE2012_CHIRP  Plot Fig. 6 or 7 of Rønne et al. (2012)
%   Usage: plotroenne2012_chirp(flag)
%
%   `plotroenne2012_chirp(waveVamp, waveVlat)` plots the output of
%   |roenne2012_chirp|_ in the style of Fig. 6 or 7 of Rønne et al. (2012).
%   Simulations are compared to data from Elberling et al. (2010).
%
%   The flag may be one of:
%
%     'twofig'   Plot the amplitude and latency in two different plot
%                windows. This is the default.
%
%     'subplot   Use subplot to position the window side-by-side
%
%     'amponly'  Plot amplitude only.
%
%     'latonly'  Plot latency only.
%
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%   See also: roenne2012_chirp, roenne2012, exp_roenne2012
%  
%   References: roenne2012modeling elberling2010evaluating zilany2007representation

% Define input flags
definput.flags.type = {'twofig','subplot','amponly','latonly'};

if 0
  definput.keyvals.waveVlat     = [6.40000000000000,6.06666666666667, ...
                      5.86666666666667;6.16666666666667,5.80000000000000, ...
                      5.50000000000000;6.10000000000000,5.73333333333333, ...
                      5.36666666666667;6.03333333333334,5.66666666666667, ...
                      5.23333333333333;5.96666666666667,5.60000000000000, ...
                      5.10000000000000;5.93333333333333,5.53333333333334, ...
                      4.56666666666667];
  definput.keyvals.waveVamp     = [0.222761519867971,0.368049605362604, ...
                      0.426295534454623;0.304245724239426, ...
                      0.613841499609571,0.662285715806521; ...
                      0.313992730990412,0.635926344773530, ...
                      0.588473779376608;0.314083457955007, ...
                      0.609949392701705,0.523100938263095; ...
                      0.308279424002208,0.569882100923284, ...
                      0.411055469911713;0.299982088057727, ...
                      0.524606196043760,0.356573917626301];
  
end;

[flags,kv]  = ltfatarghelper({}, definput,varargin);

fntz        = 12; % fontsize
coldata     = [0.7,0.7,0.7];
col         = [0,0,0];

[Delay,CElatmean,CElatstd]  = data_elberling2010('fig5');
[Delay,CEmean,CEstd]        = data_elberling2010('fig4');
Delay2                      = [Delay;Delay;Delay];  

if flags.do_subplot
  subplot(1,2,1);
end;

if flags.do_twofig
  figure;
end;

%% Plot figure 6 from Rønne et al. (2012)
if flags.do_amponly
  % Plot Elberling et al. (2010) reference data
  errorbar(Delay2',CEmean,CEstd/sqrt(20),'-', ...
           'linewidth',1.5,'color',coldata);
  
  hold all;
  
  set(gca,'fontsize',fntz);
  
  axis([-1.2 5.5 0 800]);
  xlabel('Change of delay [ms]');
  ylabel('ABR amplitude [nv]');
  text(-.7,waveVamp(1,1)*1000, '20','fontsize',fntz);
  text(-.7,waveVamp(1,2)*1000, '40','fontsize',fntz);
  text(-.7,waveVamp(1,3)*1000, '60','fontsize',fntz);
  text(-.9,waveVamp(1,3)*1000+70, 'dB nHL','fontsize',fntz);
  text(-.2,50, 'Click','fontsize',fntz);
  text(Delay(2),75, '1','fontsize',fntz);
  text(Delay(3),75, '2','fontsize',fntz);
  text(Delay(4),75, '3','fontsize',fntz);
  text(Delay(5),75, '4','fontsize',fntz);
  text(Delay(6),75, '5','fontsize',fntz);
  text(3,40, 'Chirps','fontsize',fntz);

  % Plot Simulated wave V latencies
  plot(Delay,waveVamp*1000,'-*k','linewidth',1.5);
  
  xlabel('Change of delay [ms]');
  ylabel('ABR latency [ms]');
%   axis([-1.2 6.5 0 10]);

end


if flags.do_subplot
  subplot(1,2,2);
end;

if flags.do_twofig
  figure;
end;

%% Plot figure 7 from Rønne et al. (2012)
if flags.do_latonly
    
  % Amplitude in nV
  waveVamp = waveVamp*1000; 

  set(gca,'fontsize',fntz);
  
  hold all;
  errorbar(Delay2',CElatmean,CElatstd/sqrt(20),'-','linewidth',1.5, ...
           'color',coldata)

  plot(Delay,waveVlat,'-*k','linewidth',1.5);

  % Plot text strings
  text(-.7,5.88, '60','fontsize',fntz,'color',coldata);
  text(-.7,6.59, '40','fontsize',fntz,'color',coldata);
  text(-.7,7.6, '20','fontsize',fntz,'color', coldata);
  text(-.9,8.7, 'dB nHL','fontsize',fntz,'color',coldata);
  text(5.5,waveVlat(6,1), '20','fontsize',fntz);
  text(5.5,waveVlat(6,2), '40','fontsize',fntz);
  text(5.5,waveVlat(6,3), '60','fontsize',fntz);
  text(5.3,waveVlat(6,1)+.6, 'dB nHL','fontsize',fntz);
  text(-.2,.5, 'Click','fontsize',fntz);
  text(Delay(2),1, '1','fontsize',fntz);
  text(Delay(3),1, '2','fontsize',fntz);
  text(Delay(4),1, '3','fontsize',fntz);
  text(Delay(5),1, '4','fontsize',fntz);
  text(Delay(6),1, '5','fontsize',fntz);
  text(3,.5, 'Chirps','fontsize',fntz);
  box on;
  
    axis([-2 7 0 10]);
  
end
