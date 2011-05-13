function p = langendijk(ir1,ir2,varargin)
%LANGENDIJK Localization model according to Langendijk et al. (2002)
%   Usage:    p = langendijk(ir1,ir2)
%             p = langendijk(ir1,ir2,bw,do,cp,s,bal,flow,fhigh,stim)
%
%   Input argumentss:
%     ir1     : modified impulse responses of DTFs for all positions and both ears (sorted)
%     ir2     : stored impulse responses of DTF-templates for all positions and both ears (sorted)
%
%   Output arguments:
%     p       : Predicted probability density functions for response angles with respect to target positions
%
%   LANGENDIJK(ir1,ir2,... ) results to a two dimensional matrix p.  The
%   first dimension represents all possible response positions in
%   increasing order and the second dimension all possible target
%   respectively source positions. Consequently each column describes the
%   probability density function of the response distribution for one
%   special target position. If you want to plot this matrix use
%   PLOTLANGENDIJK.
%
%   LANGENDIJK accepts the following optional parameters.
%  
%-     'bw',bw   : Bandwidth of filter bands as partial of an octave. The
%                default value is 6.
%
%-     'do',do   : Differential order. The default value is 0.
%
%-     'std'     : Use 'std' for comparison. This is the default.
%  
%-     'corr'    : Use 'corr' for comparison.
%
%-     's',s     : Standard deviation of transforming Gaussian function; default: 2
%
%-     'bal',bal : Balance of left to right channel (not included in 
%                langendijk's original comparison process); default: 1
%
%-     'flow',flow   : Start frequency of filter bank. min: 0,5kHz; default: 2kHz
%
%-     'fhigh',fhigh : End frequency of filter bank; default: 16kHz
%
%-     'stim',stim   : Applied stimulus for localization test (optional)
%
%   See also: plotlangendijk
%
%R  langendijk2002contribution
%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
  
  % FIXME: fs is hardcoded
  fs=48000;

  definput.import={'langendijkcomp'};
  definput.keyvals.bw=6;
  definput.keyvals.flow=2000;
  definput.keyvals.fhigh=16000;
  definput.keyvals.stim=[];
  
  [flags,kv]=ltfatarghelper({'bw','do','s','bal','flow','fhigh','stim'},definput,varargin);
  
  if ~isempty(kv.stim)
    ir1=halfconv( ir1,stim );
  end
  
  % model calculations
  p=zeros(size(ir2,2),size(ir1,2)); % initialisation
  
  % filter bank
  % response pdf for every target position
  for ind2=1:size(ir1,2) 
    x=averagingfb(ir1(:,ind2,:),fs,kv.flow,kv.fhigh,kv.bw);
    y=zeros(length(x),size(ir2,2),size(ir2,3)); % initialisation
    for ind=1:size(y,2) % response pdf for one target position
      y(:,ind,:)=averagingfb(ir2(:,ind,:),fs,kv.flow,kv.fhigh,kv.bw);
    end
    
    % comparison process for one target position
    p(:,ind2)=langendijkcomp(x,y,'argimport',flags,kv);
  end
  

function [ outsig ] = halfconv( ir1,stim )
% HALFCONV calculates the fast convolution with fft but without ifft
% Usage:        [ outsig ] = halfconv( ir1,stim )
% Input arguments:
%     ir1:      (modified) impulse responses of DFTs for all positions and
%               both ears
%     stim:     stimulus (time domain)
% Output argument:
%     outsig:   convolution of stimulus with DTFs in frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW
% latest update: 2010-07-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nfft = 2^nextpow2(max([size(ir1,1) length(stim)]));
  stimf=fft(stim(:),nfft);
  ir1f = fft(ir1,nfft,1);
  temp=zeros(nfft,size(ir1,2),size(ir1,3));
  for ch=1:size(ir1,3)
    for ind=1:size(ir1,2)
      temp(:,ind,ch) = ir1f(:,ind,ch).* stimf;
    end
  end
  outsig=temp;
  
  
