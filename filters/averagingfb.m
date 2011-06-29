function [ outsig ] = averagingfb( insig,fs,varargin )
%AVERAGINGFB Averaging rectangular filter bank according to Langendijk
%   Usage:        averagingfb( insig,bw,fstart,fend,fs )
%
%   Input arguments:
%      insig   : Impulse response or complex spectrum
%      fs      : Sampling frequency
%      flow    : Lowest frequency, minimum: 0.5kHz, default is 2kHz
%      fhigh   : Highest frequency, default is, default is 16kHz  
%      bw      : bandwidth, possible values 3,6,9,12, default is 6.
%
%   AVERAGINGFB(insig,fs) computes an averaging filterbank as done by
%   Langendijk et al. (2002).
%
%R  langendijk2002contribution
  
% AUTHOR : Robert Baumgartner

%   size(insig);
  
  definput.keyvals.flow=2000;
  definput.keyvals.fhigh=16000;
  definput.keyvals.bw=6;
  
  [flags,keyvals,flow,fhigh,bw]  = ltfatarghelper({'flow','fhigh','bw'},definput,varargin);
    
  % calculation
  jj=0:log(fhigh/flow)/log(2)*bw; 
  if max(imag(insig))==0
    nfft=max(4096,length(insig));
    y = abs(fft(insig,nfft,1));
  else  % input signal already in frequency domain
    y = abs(insig);
    nfft=length(insig);
  end
  
  n=round(2.^((jj)/bw)*flow/fs*nfft); %startbins
  ybw=zeros(length(jj)-1,size(insig,3)); %initialisation
  for ind=jj(1)+1:jj(end)
    nj=n(ind+1)-n(ind);
    for ch=1:size(insig,3)
      ybw(ind,ch)=sqrt(1/(nj)*sum(y(n(ind):n(ind+1)-1,ch).^2));
    end
  end
  
  outsig=20*log10(ybw);


