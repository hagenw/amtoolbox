function [waveVamp,waveVlat] = roenne2012_click(stim_level,varargin) 
%roenne2012_click  Simulate ABR respone to click
%   Usage: [waveVlat]  = roenne2012_click(stim_level).
%
%   Output parameters:
%     waveVlat   : Latency of simulated ABR wave V peak.
%     waveVamp   : Amplitude of simulated ABR wave V.
%
%   Input parameters:
%     stim_level      : Simulated levels. Default: Elberling et al. (2010)
%                       stimulus levels (20, 40, 60 dB HL) calibrated to pe
%                       SPL (+ 35.2 dB), see Rønne et al. (2012).   
%
%   `ronne2012_click(stim_level)` returns click evoked ABR wave V latencies
%   and amplitudes for a range of given stimulus levels. It simulates ABR
%   responses to click stimulus using the ABR model of Rønne et
%   al. (2012). The click stimulus is defined similar to Elberling et
%   al. (2010).
%
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%   References: roenne2012modeling elberling2010evaluating zilany2007representation
  
%definput.keyvals.stim_level     = 40:10:100;
%[flags,kv] = ltfatarghelper({},definput,varargin);

fsmod       = 200e3;        % AN model fs.
modellength = 40;           % length of modelling [ms].

% load Unitary Response
[ur,fs]=data_roenne2012;

% Output filter corresponding to recording settings.
b=fir1(200,[100/(fs/2),3000/(fs/2)]);
a=1;

%% create click stimulus
% Load click stimulus (c0) from Elberling et al. (2010).
[stim,fsstim] = data_elberling2010('stim'); 

% Define length of stimulus, uses variable modellength.
refstim = zeros(modellength/1000*fs,1);                 

% Create stimulus with chirp stimulus and concatenated zeros => combined
% length = "modellength".
refstim(1:length(stim.c0)) = stim.c0;                                       


%% Simulate ABR - loop over stimulus levels
for L = 1:length(stim_level)
  lvl = stim_level(L);
  
  % call AN model
  ANdata = zilany2007humanized(lvl, refstim, fsstim,fsmod);         
  
  % subtract 50 due to spontaneous rate. 
  ANout = ANdata'-50;                                       
  
  % Sum in time across fibers = summed activity pattern.
  ANsum = sum(ANout,2);                                     
  
  % Downsample ANsum to get fs = fs_UR = 32kHz.
  ANsum = resample(ANsum,fs,fsmod);                         
  
  % Simulated potential = UR * ANsum (* = convolved).
  simpot = filter(ur,1,ANsum);                               
  
  % apply output filter similar to the recording conditions in Elberling
  % et al. (2010).
  simpot = filtfilt(b,a,simpot);                             
  
  % Find max peak value (wave V).
  maxpeak = max(simpot);                                      
  
  % Find corresponding time of max peak value (latency of wave V).
  waveVlat(L)= find(simpot == max(simpot));              
  
  % find minimum in the interval from "max peak" to 6.7 ms later.
  minpeak = min(simpot(find(simpot == maxpeak):find(simpot == maxpeak)+100)); 
  
  % Calculate wave V amplitude, as the difference between the peak and
  % the following dip.
  waveVamp(L) = (maxpeak-minpeak);                           
  
end

% Subtract 15 ms as click stimulus peaks 15 ms into the time series
waveVlat = waveVlat/fs*1000-15;