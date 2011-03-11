%DEMO_ABSOLUTETHRESHOLD  Plot of the different absolute thresholds
%
%  This demos generates a simple figure that shows the behaviour of
%  the different absolute thresholds of hearing.
%
%  FIGURE 1  Thresholds of hearing by standard
%
%      This demos generates a simple figure that shows the behaviour of
%      the different absolute thresholds of hearing.    
%
%  See also:  absolutethreshold

disp(['Type "help demo_absolutethreshold" to see a description of how this ', ...
      'demo works.']);

% Set the limits
flow=125;
fhigh=8000;
plotpoints=50;

xrange=linspace(flow,fhigh,plotpoints);


figure(1)

types   = {'iso226_2003','map','er3a','er2a'};
symbols = {'k-'         ,'ro' ,'gx'  ,'b+'  };

hold on;
for ii=1:numel(types)
  semiaudplot(xrange,absolutethreshold(xrange,types{ii}),'opts',{symbols{ii}});
end;
hold off;
legend(types{:},'Location','SouthEast');
xlabel('Frequency (Hz)');
ylabel('Absolte threshold (dB SPL)');

