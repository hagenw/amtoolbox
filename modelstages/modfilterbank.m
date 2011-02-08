function [outsig,mfc] = modfilterbank(insig,fs,fc,varargin)
%MODFILTERBANK  Modulation filter bank
%   Usage: [outsig, mfc] = modfilterbank(insig,fs,fc);
%

% fs    = sampling rate in Hz,
%         should be greater than 8000 Hz to avoid aliasing errors.
% 
% coef	= center frequencies of the modulation filters.
% [out1,out2, ...,outn] = each column of martrix out contains the output of
%         a single modulation filter.
%
% copyright (c) 1999 Stephan Ewert and Torsten Dau, Universitaet Oldenburg

nfreqchannels=length(fc);

Q = 2;
bw = 5;
ex=(1+1/(2*Q))/(1-1/(2*Q));

outsig=cell(nfreqchannels,1);

startmf = 5;

% second order modulation Butterworth lowpass filter with a cut-off frequency of 2.5
% Hz.
[b_lowpass,a_lowpass] = solp(2*pi*2.5/fs,1/sqrt(2));

% first order modulation Butterworth lowpass filter with a cut-off
% frequency of 150 Hz. This is to remove all modulation frequencies
% above 150 Hz.
[b_highest,a_highest] = butter(1,150/(fs/2));

% Set the highest modulation frequency as proportion of the corresponding
% center frequency.
umf = min(fc.*0.25, 1000);  

for freqchannel=1:nfreqchannels

  % Cut away highest modulation frequencies
  outtmp = filter(b_highest,a_highest,insig(:,freqchannel));

  if umf(freqchannel)==0
    % ----------- only lowpass ---------------------
    outsigblock = filter(b_lowpass,a_lowpass,outtmp);
    mfc = 0;

  else                
    tmp = fix((min(umf(freqchannel),10) - startmf)/bw);
    tmp = 0:tmp;
    mfc = startmf + 5*tmp;
    tmp2 = (mfc(end)+bw/2)/(1-1/(2*Q));
    tmp = fix(log(umf(freqchannel)/tmp2)/log(ex));
    tmp = 0:tmp;
    tmp = ex.^tmp;
    mfc=[mfc tmp2*tmp];

    % --------- lowpass and modulation filter(s) ---
    outsigblock = zeros(length(insig),length(mfc)+1);
    outsigblock(:,1) = filter(b_lowpass,a_lowpass,outtmp);
    mfc = [0 mfc];

    for nmfc=2:length(mfc)
      w0 = 2*pi*mfc(nmfc)/fs;
      if mfc(nmfc) < 10
   	[b3,a3] = efilt(w0,2*pi*bw/fs);
      else
        [b3,a3] = efilt(w0,w0/Q);
      end
      
      outsigblock(:,nmfc) = 2*filter(b3,a3,outtmp);
    end
    
  end
  
  %% ------------ post-processing --------------------
  
  for nmfc=1:length(mfc) % v2 MJ 17. oct 2006
    if mfc(nmfc) <= 10
      outsigblock(:,nmfc) = 1*real(outsigblock(:,nmfc));
    else
      outsigblock(:,nmfc) = 1/sqrt(2)*abs(outsigblock(:,nmfc));
    end
  end
  

  outsig{freqchannel}=outsigblock;
end;

%% ------------ subfunctions ------------------------


% complex frequency shifted first order lowpass
function [b,a] = efilt(w0,bw);

e0 = exp(-bw/2);

b = 1 - e0;
a = [1, -e0*exp(1i*w0)];

% second order Butterworth lowpass filter
function [b,a] = solp(w0,Q)

W0 = tan(w0/2);

b = [1; 2; 1];
a = [1 + 1/(Q*W0) + 1/W0^2; -2/W0^2 + 2; 1 - 1/(Q*W0) + 1/W0^2];

b = b/a(1);
a = a/a(1);

% first order lowpass filter (from mfb2.m - MJ)
function [b,a] = folp(w0);

W0 = tan(w0/2);

b = [W0, W0]/(1 + W0);
a = [1,(W0 - 1)/(1 + W0)];

% end of mfbtd.m







