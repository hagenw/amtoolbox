function definput=arg_baumgartner2014(definput)
  
definput.flags.regularization = {'regular','noregular'};
definput.flags.motoricresponsescatter = {'mrs','nomrs'};
definput.flags.settings = {'notprint','print'};

% CP-Falgs:
definput.flags.cp={'fwstd','std','xcorr'};

definput.keyvals.fs=48000;      % sampling rate in Hz
definput.keyvals.S=0.5;         % listener-specific sensitivity
definput.keyvals.lat=0;         % lateral angle in deg
definput.keyvals.stim=[];
definput.keyvals.fsstim=[];
definput.keyvals.space=1;       % No. of ERBs (Cams) 
definput.keyvals.flow=700;      % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.do=1;          % differential order
definput.keyvals.bwcoef=13;     % binaural weighting coefficient in deg
definput.keyvals.polsamp=[-30:5:70 80 100 110:5:210];  % polar sampling (for regular flag)
definput.keyvals.rangsamp=5;    % equi-polar sampling of response angles
definput.keyvals.mrsmsp=17;     % response scatter (in deg) in the median plane induced by non-auditory processes
definput.keyvals.gamma=6;       % degree of selectivity in 1/dB
  
end