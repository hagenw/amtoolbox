function [out, CenterFreq, MF_CFs] = ref_dau1997preproc(insig, fs)
%
%  This is a trimmed down version of CASP_Preproc by Morten Løve Jepsen
%  to just reproduce the dau1997 model. All values are hardwired.

CFlow=80;
CFhigh=8000;
baseF=1000;
  
inLen = length(insig);

% find the center frequencies used in the filterbank, 1 ERB spacing
[NrFBChannels, CenterFreq] = getFBCenterERBs(CFlow, CFhigh, baseF);

%  GT coeffs, SE
[GT_b(1,:,:), GT_a(1,:,:)]=getGFBFilterCoefs(NrFBChannels, CenterFreq, 1, fs);

% lowest and highest CFs of the MFB as function of CF
MFlow = CenterFreq .* 0;                        % set lowest mf as constant value
MFhigh = min(CenterFreq .* 0.25, 1000);         % set highest mf as proportion of CF
[MF_CFs,out] = mfbtd(1,min(MFlow),max(MFhigh),1,fs); % to find the number of MF's

NrMFChannels = size(MF_CFs,2);                  % maximum number of modulation filters

out = zeros(inLen,NrFBChannels,NrMFChannels); % define output array

for ChannelNr = 1:NrFBChannels
    
  % Gammatone filterbank
  current_GTb = squeeze(GT_b(1,ChannelNr,:));
  current_GTa = squeeze(GT_a(1,ChannelNr,:));
  y = 2*real(filter(current_GTb,current_GTa,insig));
  
  % 'haircell' envelope extraction
  y = max( y, 0 );
  [LP1000_b, LP1000_a] = butter(2, 1000*2/fs);
  y = filter(LP1000_b,LP1000_a, y);   % 2nd order butterworth
  y = nlal_lim(y, fs,10);             % AD loops, 15 is overshoot limit factor
  
  % Modulation filterbank
  [infpar,y] = mfbtd(y,MFlow(ChannelNr),MFhigh(ChannelNr),1,fs);	% MFB incl 150 LP
  y = mfbtdpp(y,infpar,fs);
  
  % Fill 'y' into output array 
  out(:,ChannelNr,1:length(infpar)) = y;
            
end


