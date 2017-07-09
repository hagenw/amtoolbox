function pmv = langendijk2002(targets,template,varargin)
%LANGENDIJK2002 Localization model according to Langendijk et al. (2002)
%   Usage:    pmv = langendijk2002(targets,template)
%             pmv = langendijk2002(targets,template,fs,bw,s,do,flow,fhigh)
%
%   Input parameters:
%     targets  : head-related impulse responses (HRIRs) of target sounds 
%                (sorted acc. ascending polar angle)
%     template : HRIRs of template
%
%   Output parameters:
%     pmv     : Predicted probability mass vectors (PMVs) of polar response
%               angles as a function of the polar target angle.
%
%   `langendijk2002(targets,template,... )` results to a two dimensional matrix p.  The
%   first dimension represents all possible response positions in
%   increasing order and the second dimension all possible target
%   respectively source positions. Consequently each column represents the
%   predicted probability mass vector (PMV) of the polar response angle 
%   distribution for one special target position. If you want to plot this 
%   prediction matrix use |plot_langendijk2002|.
%
%   `langendijk2002` accepts the following optional parameters.
%
%     'fs',fs        Sampling rate of the head-related impulse responses.
%  
%     'bw',bw        Bandwidth of filter bands as partial of an octave. The
%                    default value is 6.
%
%     'do',do        Differential order. The default value is 0.
%
%     's',s          Standard deviation of transforming Gaussian
%                    function; default value is 2.
%
%     'flow',flow    Lower cutoff frequency of filter bank. min: 0,5kHz; default: 2kHz
%
%     'fhigh',fhigh  Upper cutoff frequency of filter bank; default: 16kHz
%
%   `langendijk2002` accepts the following flags.
%
%     'std'          Apply Gaussian transformed standard deviation of 
%                    inter-spectral differences for comparison process. 
%                    This is the default.
%  
%     'xcorr'        Apply crosscorrelation for comparison process.
%
%   See also: plot_langendijk2002
%
%   References: langendijk2002contribution

% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
  
  
  definput.import={'langendijk2002_comp'};
  definput.keyvals.bw=6;
  definput.keyvals.flow=2000;
  definput.keyvals.fhigh=16000;
  definput.keyvals.stim=[];
  definput.keyvals.fs=48000;
  
  [flags,kv]=ltfatarghelper({'fs','bw','s','do','flow','fhigh'},definput,varargin);
  
  % Stimulus (not considered in original model)
  if not(isempty(kv.stim))
    tmp = convolve(kv.stim,targets);
    targets = reshape(tmp,[size(tmp,1),size(targets,2),size(targets,3)]);
  end
  
  % Filter bank
  x = cqdft(targets,kv.fs,kv.flow,kv.fhigh,kv.bw);
  y = cqdft(template,kv.fs,kv.flow,kv.fhigh,kv.bw);
  
  % Comparison process
  si=zeros(size(template,2),size(targets,2),size(template,3)); % initialisation
  for ii=1:size(targets,2)
      si(:,ii,:) = langendijk2002_comp(x(:,ii,:),y,'argimport',flags,kv);
  end
  
  % Binaural average
  si = mean(si,3);
  
  % Normalization to PMV
  pmv = si ./ repmat(sum(si),size(si,1),1);
  
end