%DEMO_ABSOLUTETHRESHOLD  Standards for absolute threshold of hearing
%
%  This demos generates a simple figure that shows the behaviour of
%  the different standards for absolute thresholds of hearing.
%
%  FIGURE 1  Thresholds of hearing by standard
%
%      This demos generates a simple figure that shows the behaviour of
%      the different absolute thresholds of hearing.    
%
%  FIGURE 2 High frequencies
%
%      Absolute thesholds for the ER2A and the Sennheiser HDA-200
%      earphones are provided up to 16 kHz.
%
%  See also:  absolutethreshold

disp(['Type "help demo_absolutethreshold" to see a description of how this ', ...
      'demo works.']);


figure(1);

flow=125;
fhigh=8000;
plotpoints=30;

xrange=linspace(flow,fhigh,plotpoints);

types   = {'iso226_2003','map','er3a','er2a','hda200'};
legends = {'iso226\_2003','map','er3a','er2a','hda200'};
symbols = {'k-'         ,'ro' ,'gx'  ,'b+','y*'  };

hold on;
for ii=1:numel(types)
  semiaudplot(xrange,absolutethreshold(xrange,types{ii}),'opts',{symbols{ii}});
end;
hold off;
legend(legends{:},'Location','SouthEast');
xlabel('Frequency (Hz)');
ylabel('Absolte threshold (dB SPL)');


figure(2);

flow=125;
fhigh=16000;
plotpoints=30;

xrange=linspace(flow,fhigh,plotpoints);

types   = {'er2a','hda200'};
legends = {'er2a','hda200'};
symbols = {'b+','y*'  };

hold on;
for ii=1:numel(types)
  semiaudplot(xrange,absolutethreshold(xrange,types{ii}),'opts',{symbols{ii}});
end;
hold off;
legend(legends{:},'Location','SouthEast');
xlabel('Frequency (Hz)');
ylabel('Absolte threshold (dB SPL)');


%OLDFORMAT
