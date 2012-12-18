function [out, CenterFreq, MF_CFs] = ref_jepsen2008preproc(insig, fs)
%
%  This is a trimmed down version of CASP_Preproc by Morten Løve Jepsen
%  to just reproduce the jepsen2008. All values are hardwired.
  
CFlow=80;
CFhigh=8000;
baseF=1000;

%
% usage [out, CenterFreq, MF_CFs] = CASPHI_Preproc(in, fs, CFlow, CFhigh,
% baseF, BMtype, MFtype, subj)
%   in = input acoustic signal, 
%   fs = sampling rate, 
%   CFlow   = lowest filter center frequency, 
%   CFhigh  = highest CF, 
%   baseF   = base frequency for filterbank,
%   BMtype: 'gt' = gammatone or 
%           'drnl' = dual-resonance NL
%   MFtype: 'lp' = modulation lowpass filter as in Dau et al 1996 or
%           'mfb' = for the modulation filterbank modified in Jepsen 07
%   subj:   'NH' = normal hearing
%           'HIx' = hearing impaired with lost compression and 20 dB flat
%           IHC loss
%
% 2. feb 2009, Morten Løve Jepsen

%function [out, CenterFreq, MF_CFs] = CASP_Preproc(in, fs, CFlow, CFhigh, baseF, BMtype, MFtype)

inLen = length(insig);

% find the center frequencies used in the filterbank, 1 ERB spacing
[NrFBChannels, CenterFreq] = getFBCenterERBs(CFlow, CFhigh, baseF);


%  Filter with outer and middle ear TF
x_stapes = OuterMiddleFilter(insig); %note that this one only applyes to 44.1 kHz sampling
load minlim.mat

% lowest and highest CFs of the MFB as function of CF
MFlow = CenterFreq .* 0;                        % set lowest mf as constant value
MFhigh = min(CenterFreq .* 0.25, 1000);         % set highest mf as proportion of CF
[MF_CFs,out] = mfbtd(1,min(MFlow),max(MFhigh),1,fs); % to find the number of MF's

NrMFChannels = size(MF_CFs,2);                  % maximum number of modulation filters

out = zeros(inLen,NrFBChannels,NrMFChannels); % define output array

for ChannelNr = 1:NrFBChannels
    
  % DRNL filterbank
  % x_stapes = OuterMiddleFilter(insig,fs);
  y = drnl(x_stapes',CenterFreq(ChannelNr),fs)';
  
  % 'haircell' envelope extraction
  y = max( y, 0 );
  [LP1000_b, LP1000_a] = butter(2, 1000*2/fs);
  y = filter(LP1000_b,LP1000_a, y);   % 2nd order butterworth
                                      %%%%%%%
  % Expansion and compensate for middle ear TF
  y = y*10^(50/20);                   % linear gain to fit ADloop
  y = y.^2;                           %expansion
           
  % non-linear adaptation loops
  y = max(y,8e-5);
  y = nlal_lim(y, fs,15);             % AD loops, 15 is overshoot limit factor

  % Modulation filterbank
  [infpar,y] = mfbtd(y,MFlow(ChannelNr),MFhigh(ChannelNr),1,fs);	% MFB incl 150 LP
  y = mfbtdpp(y,infpar,fs);
  
  % Fill 'y' into output array 
  out(:,ChannelNr,1:length(infpar)) = y;
            
end


