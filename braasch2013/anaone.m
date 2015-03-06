function [xLgS,xLgL,lag,B,midFreq]=anaone(x,Fs);

% function [xb,x2a,x2b,lag,B,midFreq]=anaone(x,Fs);

% Monaural preprocessing (Section II.E) &
% Lag removal (Section II.B) for acmod
%
% (c) 2013, Rensselaer Polytechnic Institute
% contact: jonasbraasch@gmail.com
%
% INPUT PARAMETERS:
% y        = monaural left or right input signal
% Fs       = sampling frequency (tested for 48 kHz)
% OUTPUT PARAMETERS:
% xLgS    = signal with removed lag assuming lag smaller than lead
% xLgL    = signal with removed lag assuming lag larger than lead
% lag     = lag delay in tabs
% B       = relative lag amplitude according to Eq. 2
% midFreq = Center frequencies for auditory bands


Fms=Fs./1000; % Sampling frequency based on milliseconds
lags=round(10.*Fms); % number of lags for autocorrelation process 
MainPeakWidth=0.75.*Fms; % max main peak width to determine minimum lag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Monaural Preprocessing (Section II.E)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auditory filterbank from AMT toolbox
[xb,midFreq]=auditoryfilterbank(x,Fs,'flow',100,'fhigh',1400);

% Halfwave rectification and then correction
for n=1:length(midFreq) % loop over all auditory bands
    index=find(xb(:,n)<0); % find values below zero
    xb(index,n)=0; % set these values to zero (halfwave rectification)
    
    % reconstruct halfwave based on center frequency of auditory band
    % see Eq. 6 and Fig. 8
    xb(:,n)=reconHW(xb(:,n),Fs,midFreq(n)); 
end % of for    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Lag Removal (Section II.B)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine autocorrelation function for lag detection, see Eq. A2
for n=1:length(midFreq) % loop over all auditory bands
    xba(n,:)=xcorr(xb(:,n)./(max([midFreq(n) 100])).^0.5,xb(:,n)./(max([midFreq(n) 100]))^0.5,lags)';
end

% BEGIN amplitude normalization across auditory bands
% see Eqs. 7 and 8
% necessary to simulate psychoacoustical results for narrowband signals
% can be commented out for broadband signals

maxum=20*log10(max(xb));
maxum=maxum-20*log10(max(max(xb)));
index=find(maxum>-12);
bandwidth=length(index);

for n=1:length(midFreq) % loop over all auditory bands
        if bandwidth<=5; 
            maxall=max(max(xba));
            maxi=max(xba(n,:)); % maximum in current auditory band
            % find bands with energy greater than -60dB from maximum 
            if 20*log10(maxi./maxall)>-60 
                xba(n,:)=xba(n,:)./maxi; % normalize to one
            else
                xba(n,:)=xba(n,:).*0; % set to zero
            end % of 
        end % of if 
    
end % of for
% END amplitude normalization across auditory bands

% Now determine lag parameters from frequency-based 
% autocorrelation functions xba 
xa=sum(xba); % integrated autocorrelation function with normalized functions
Tmaxi=max(xa); % find main maximum of autocorrelation function

% autocorrelation function right of the main peak to exclude main peak
xaSide=xa(lags+MainPeakWidth+1:lags*2+1); 
[maxi,index]=max(xaSide); % find side peak based on maximum
lag=index(1)+MainPeakWidth-1; % delay between lead and lag, Sec. II.A, Fig. 2
B=PeakRatio(maxi/Tmaxi); % amplitude ratio between lead and lag, Eq. 2

% Limit B to 0.9, see Sec. II.C, p. 423
if B>0.9 
    B=0.9;
end % of if 

% Took precise calculation for very small lags (<1ms at contralateral side out)
% Will need to do this after level-calibrating new filterbank. Will only
% affect 1-ms ISI condition or smaller ISIs. 

% Now, we will remove the lag twice: (1) assuming that the lag amplitude is
% smaller than the lead amplitude, (2) assuming that the lag amplitude is
% larger than the lead amplitude
% We, will later detemine binaurally in function acmod.m which was the
% correct assumption using the LLAR modes.

for n=1:length(midFreq) % loop over all auditory bands
    xLgS(n,:)=de_conv(xb(:,n),lag,B,Fs)'; % (1) deconvolution with smaller lag amp
    xLgL(n,:)=rev_conv(xb(:,n),lag,B,Fs)'; % (2) deconvolution with larger lag amp
end % of for