function [ITD,ILD,lagL,lagR,BL,BR,LLARmode]=acmod(y,Fs);

% function [ITD,ILD,lagL,lagR,BL,BR,LLARmode,info]=acmod(y,Fs);
%
% Implementation of a precedence effect model to simulate localization 
% dominance using an adaptive, stimulus parameter-based inhibition process
% 
% J. Braasch (2013) J. Acoust. Soc. Am. 134(1): 420-435.
% For section references see paper
%
% (c) 2013, Rensselaer Polytechnic Institute
% contact: jonasbraasch@gmail.com
%
% INPUT PARAMETERS:
% y        = binaural input signal, y=y(samples,stereo channels)
% Fs       = sampling frequency (tested for 48 kHz)
% OUTPUT PARAMETERS:
% ITD      = Interaural Time Difference collapsed over time and frequency
% ILD      = Interaural Level Difference collapsed over time and frequency
% lagL     = lag delay for left lag [ms]
% lagR     = lag delay for right lag [ms]
% BL       = left-lag amplitude, relative according to Eq. A5    
% BR       = right-lag amplitude, relative according to Eq. A5    
% LLARmode = Lag-to-Lead Amplitude Ratio mode, see Table 1
%
% dependent functions:
% db2amp.m    =
% anaone.m    =
% de_conv.m   =
% derive.m    =
% ncorr.m     =
% rev_conv.m  =
% PeakRatio.m =
% reconHW.m   =
%
% Funtions used from AMT toolbox: auditoryfilterbank.m

Fms=Fs./1000; % sampling frequency based on milliseconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monaural preprocessing (Section II.E) ^
% Lag removal (Section II.B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yL_LgS,yL_LgL,lagL,BL,midFreq]=anaone(y(:,1),Fs); % left channel
[yR_LgS,yR_LgL,lagR,BR,midFreq]=anaone(y(:,2),Fs); % right channel
% yL_LgS=left signal with removed lag assuming lag smaller than lead
% yL_LgL=left signal with removed lag assuming lag larger than lead
% yR_LgS=right signal with removed lag assuming lag smaller than lead
% yR_LgL=right signal with removed lag assuming lag larger than lead

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine correct LLAR mode
% Section II.D, Table 1 & Fig. 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cross-correlation to determine best combination
n1=max(ncorr(sum(yL_LgS),sum(yR_LgS),Fms)); % calc ICC ofr LLAR mode 1
n2=max(ncorr(sum(yL_LgL),sum(yR_LgL),Fms)); % calc ICC ofr LLAR mode 2
n3=max(ncorr(sum(yL_LgS),sum(yR_LgL),Fms)); % calc ICC ofr LLAR mode 3
n4=max(ncorr(sum(yL_LgL),sum(yR_LgS),Fms)); % calc ICC ofr LLAR mode 4
% Pick LLARmode based on highest coherence
[maxi,LLARmode]=max([n1 n2 n3 n4]);

% asign left and right signals x1/x2 form best LLAR mode
switch LLARmode
    case 1
        x1=yL_LgS;
        x2=yR_LgS;
    case 2
        x1=yL_LgL;
        x2=yR_LgL;
    case 3
        x1=yL_LgS;
        x2=yR_LgL;
    case 4
        x1=yL_LgL;
        x2=yR_LgS;
end % switch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of binaural cues in individual frequency bands
% Section III.A, see Box "CE” in Fig. 9
%
% Note: In this model version we do NOT select LLAR Mode 1 
% automatically if all four correlation values
% differ by less than 0.1 (p. 428)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

windowlength=50.*Fms; % windowlength for temporal steps
steps=ceil(length(x1(1,:))/windowlength); % determine number of steps
h=triang(windowlength); % window for overlap-add method

for n=1:length(midFreq) % loop over all frequency bands
    xL=[x1(n,:)'; zeros(windowlength*2,1)]; % add zeros to analyze whole file
    xR=[x2(n,:)'; zeros(windowlength*2,1)]; % add zeros to analyze whole file
    for m=0:steps*2-2 % loop over all time steps
        %window out left and right signals
        yL=xL(1+m*round(windowlength/2):windowlength+m*round(windowlength/2)).*h;
        yR=xR(1+m*round(windowlength/2):windowlength+m*round(windowlength/2)).*h;
        ICC(m+1,:)=xcorr(yL,yR,Fms)'; % interaural cross-correlation
        eL=mean(sqrt(yL.^2)); % rms energy in left channel
        eR=mean(sqrt(yR.^2)); % rms energy in right channel
        if eL>0 & eR>0; % if signal in both channels exist
            ILD_f(m+1)=20*log10(eL./eR); % ILD within frequency band
            En(m+1)=eL.*eR; % Amplitude
        else 
            ILD_f(m+1)=0;
            En(m+1)=0;
        end % of of                    
    end % of for 

    % Frequency weighting according to Stern et al. (1988)
    % frequency weighting filter coefficients: b1, b2, b3
    b1=-0.0938272;
    b2=0.000112586;
    b3=-0.0000000399154;
    f=midFreq(n);
    if f<1200
        w=10.^((-b1.*f-b2.*f.^2-b3.*f.^3)./10)./276;
    else
        w=10.^((-b1.*1200-b2.*1200.^2-b3.*1200.^3)./10)./276;
    end % of if 
    ICCintT(n,:)=sum(ICC).*w; % integrate ICC over time with Freq weighting
    ILDint(n)=sum(ILD_f.*En)./sum(En); % integrated ILD, energy weighted
    Eint(n)=sum(En); % integrated Energy
end % of for 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of binaural cues across frequency bands
% See Section III.A.1. Decision device
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ICCintF=sum(ICCintT); % Detemine ICC collapsed over all frequency bands

% Normalize signal for centroid calculation
ICCintF=ICCintF-min(ICCintF);
index=find(ICCintF<max(ICCintF)./2);
ICCintF(index)=0;

ITD=sum((-Fms:Fms).*ICCintF./sum(ICCintF))./Fms; % ITD estimation based on centroid

% ILD calculation amplitude weighted over frequency bands
ILD=sum(ILDint.*Eint)./sum(Eint);

% convert lag delays to ms:
lagL=lagL./Fms;
lagR=lagR./Fms;

