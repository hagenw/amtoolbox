function [waveVamp, waveVlat]  = roenne2012(stim,fsstim,stim_level,varargin)
%ROENNE2012 Simulates an ABR to any given stimulus
%   Usage: [waveVamp, waveVlat]  = roenne2012(flag)
%
%   Output parameters:
%     waveVamp   : Amplitude of simulated ABR wave V.
%     waveVlat   : Latency of simulated ABR wave V peak.
%
%   `roenne2012(stim,fsstim,stim_level)` returns simulated ABR wave V
%   latency and amplitude. The stimulus *stim* must be defined in pascals
%   and calibrated so a pure tone stimulus has an RMS value of 1. Transient
%   stimuli (which this model is designed to simulate) has to be calibrated
%   in peSPL acoustically. This is **not** the same as "just" having a
%   numerical peak to peak value of the same level as the pure tone. For
%   calibrated click, chirps and tone bursts, see |roenne2012_click|,
%   |roenne2012_tonebursts| and |roenne2012_chirp|.
%
%   The parameter *fsstim* gives the sampling frequency of the input
%   stimulus, and *stim_level* the level. As input is calibrated to an
%   RMS-value of 1, a stimulus level in (pe)SPL has to be set.
%
%   The flag may be one of:
%
%     'plot'            Plot the output. See |plot_roenne2012|.
%  
%     'noplot'          Do not plot. This is the default.
%
%     'fsmod',fsmod     Auditory nerve model sampling frequency.
%                       Default value is 200000.
%      
%     'flow',flow       Auditory nerve model lowest center frequency.
%                       Default value is 100 Hz.
%
%     'fhigh',fhigh     Auditory nerve model highest center frequency.
%                       Default value is 16000 Hz.
%
%     'min_modellength',mn 
%                       Minimum length of modelling measured in ms.
%                       Default value is 40.
%
%   Examples:
%   ---------
%
%   Simulates a click evoked ABR (c0 of the loaded file is a click). Note
%   that the click loaded in this example starts after 15ms. The simulated
%   wave V latency is thus also 15 ms "late" :::
%
%     stim=data_elberling2010('stim'); 
%     roenne2012(stim.c0,30e3,60,'plot')
%
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%   References: roenne2012modeling elberling2010evaluating zilany2007representation

% Define input flags
definput.flags.plot     = {'plot','noplot'};
definput.keyvals.fsmod=200000;
definput.keyvals.flow = 100;
definput.keyvals.fhigh = 16000;
definput.keyvals.min_modellength=40;
[flags,kv]      = ltfatarghelper({},definput,varargin);

%% Init
[ur,fs] = data_roenne2012;

% Assure minimum model length of 40ms
if length(stim)/fsstim < kv.min_modellength/1000                               
    stim_temp = zeros(1, fsstim*kv.min_modellength/1000);
    stim_temp(1:length(stim)) = stim;
    stim = stim_temp;
end

%% ABR model
% call AN model, note that lots of extra outputs are possible
[ANout,vFreq] = zilany2007(stim_level, stim, fsstim, kv.fsmod, 'flow',kv.flow, 'fhigh',kv.fhigh);   

% subtract 50 due to spontaneous rate
ANout = ANout'-50;                                            

% Sum in time across fibers, summed activity pattern
ANsum1 = sum(ANout,2);                 

% Downsample ANsum to get fs = fs_UR = 32kHz
ANsum = resample(ANsum1,fs,kv.fsmod); 

% Simulated potential = UR * ANsum (* = convolution)
simpot = filter(ur,1,ANsum);        

% Find max peak value (wave V)
maxpeak = max(simpot);                                          

% Find corresponding time of max peak value (latency of wave V). The unit
% is [ms]. 
waveVlat = find(simpot == maxpeak)/fs*1000; 

% find minimum in the interval from "max peak" to 6.7 ms later
minpeak = min(simpot(find(simpot == max(simpot)):...
                     find(simpot == max(simpot))+200)); 

% Calculate wave V amplitude, as the difference between the peak and the
% dip, in [\mu p] (micro pascals).
waveVamp = (maxpeak-minpeak);                               

if flags.do_plot
  plot_roenne2012(stim_level,waveVamp, waveVlat, simpot, ANout, 'flow',kv.flow, 'fhigh', kv.fhigh);
end
