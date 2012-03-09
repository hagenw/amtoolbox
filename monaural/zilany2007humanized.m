function [ANdata,vFreq] = zilany2007humanized(stim_level,stim,fswav,fsmod,varargin)
% ZILANY2007HUMANIZED  Humanized auditory nerve model
%   Usage: [ANdata,vFreq] = zilany2007humanized(lvl,stim,fswav,fsmod);
%
%   Input parameters:
%     fswav       : Sampling frequency of stimulus
%     stim_level  : Level of stimulus in peSPL
%     stim        : Pressure waveform of stimulus (timeseries)
%     fsmod       : Model sampling frequency (often 200kHz)
%
%   Output parameters:
%     ANdata     : AN exicitation in 500 different AN fibers spaced equally
%                   on the BM
%     vFreq      : Frequency vector containing the 500 center frequencies
%
%   `zilany2007humanized(lvl, stim, fswav, fsmod)` returns simulations from
%   Rønne et al. (2012). It calls the mex'ed C code containing the humanized
%   version of Zilany et al. (2007)'s AN model. The humanization is
%   described in Rønne et al. (2012). The AN model is called 500 times to
%   simulate 500 fibers tuned to different center frequencies.
%
%   Please cite Rønne et al (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
%   This function takes the following optional parameters:
%
%     'flow',flow     Lowest centre frequency. Default value is 100.
%
%     'fhigh',fhigh'  Highest centre frequency. Default value is 16000.
%
%     'nfibers',nf    Number of fibers between lowest and highest
%                     frequency. The fibers will be equidistantly spaced
%                     on the basilar membrane. Default value is 500.
%  
%   References: roenne2012modeling zilany2007representation

if nargin<4
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Define input flags
definput.keyvals.flow    = 100;
definput.keyvals.fhigh   = 100;
definput.keyvals.nfibers = 16000;
[flags,kv]  = ltfatarghelper({'flow','fhigh','nfibers'},definput,varargin);

fs      = fsmod;                                    % model fs
stim    = resample(stim,fs,fswav);                  % stim fs = mod fs
stim    = 20e-6*10^(stim_level/20)*stim;            % Calibrate level

% stim must be a row vecto
if size(stim,2) == 1
    stim = stim';
end

% number of fibres between cflo and cfhi
NFibers = 500;                                

% location of lowest and highest centre frequency
xlo     = (1.0/0.06)*log10((kv.flow/165.4)+0.88);
xhi     = (1.0/0.06)*log10((kv.fhigh/165.4)+0.88);	

% equal spaced distances on the BM
vX      = linspace(xlo,xhi,kv.nfibers);

% and the resulting frequency vector
vFreq   = 165.4*(10.^(0.06*vX)-0.88); 

% resolution in the time domain
tdres   = 1/fs;                

% spontaneous rate in sp/sec
spont   = 50;                  

% cohc is the ohc scaling factor: 1 is normal OHC function; 0 is complete
% OHC dysfunction
cohc    = 1;                                        

% cihc is the ihc scaling factor: 1 is normal IHC function; 0 is complete
% IHC dysfunction
cihc    = 1;                                        

% Call AN model - loop over the 500 fibers tuned to different CFs
for jj = 1:kv.nfibers
  
  % Call AN model (mex'ed C model)
  [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth500k] ...
      = comp_zilany2007humanized(stim,...
                                 vFreq(jj),...
                                 1,...
                                 tdres,...
                                 length(stim)/fs,cohc, ...
                                 cihc,...
                                 spont); 
  
  % Use the output of the synapse stage of the AN model. 
  ANdata(jj,:)           = synout;             
end