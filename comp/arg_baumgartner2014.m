function definput=arg_baumgartner2014(definput)
  
definput.flags.fbank = {'gammatone','cqdft','drnl','zilany2007humanized','zilany5'};   % disp('zilany2007humanized used!')
definput.flags.headphonefilter = {'','headphone'};
definput.flags.middleearfilter = {'','middleear'};
definput.flags.ihc = {'noihc','ihc'};    
definput.flags.Ifw = {'nointensityweighting','intensityweighting'};
definput.flags.regularization = {'regular','noregular'};
definput.flags.motoricresponsescatter = {'mrs','nomrs'};
definput.flags.sensitivitymapping = {'sigmoidmapping','normstdmapping'};
definput.flags.settings = {'notprint','print'};

% CP-Falgs:
definput.flags.cp={'fwstd','std','xcorr'};

definput.keyvals.fs=48000;      % Hz
definput.keyvals.S=0.5;         % listener-specific sensitivity parameter
definput.keyvals.lat=0;         % deg
definput.keyvals.stim=[];
definput.keyvals.fsstim=[];
definput.keyvals.space=1;       % No. of ERBs (Cams) 
definput.keyvals.do=1;
definput.keyvals.flow=700;      % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.lvltar = 50; 	% dBSPL
definput.keyvals.lvltem = 60;  	% dBSPL
definput.keyvals.SL = [];       % db/ERB; spectral density of target sound re absolut detection threshold
definput.keyvals.conGain=13;    % steepness in degrees of binaural weighting function
definput.keyvals.polsamp=[-30:5:70 80 100 110:5:210];  % polar sampling (for regular)
definput.keyvals.mrsmsp=19;     % degrees
definput.keyvals.gamma=6;       % slope of psychometric function
  
end