function exp_georganti2013(varargin)
%EXP_GEORGANTI2013 Figures from Georganti et al. (2013)
%   Usage: exp_georganti2013(flag);
% 
%   `exp_georganti2013(flags)` reproduces figures from Georganti et al. (2013).
%   
%   The following flags can be specified;
%
%     'fig9'    distance estimation for a small room
%
%     'fig10'    distance estimation for a large room
%
%   The figures present the computed distance-dependent feature
%   BSMD STD (Binaural Spectral Magnitude Difference Standard Deviation)
%   for audio signals measured at different/receiver distances.
%
%   The BSMD STD may be derived from any dual-channel signal
%   (binaural/stereo recordings).
%
%   The BSMD STD feature is related to the standard deviation of the
%   magnitude spectrum of room impulse response, which is known to depend on 
%   the source/receiver distance. See [2] for more information.
%
%
%   REFERENCES:
%
%   [1] E. Georganti, T. May, S. van de Par, and J. Mourjopoulos. "Extracting
%   sound-source-distance information from binaural signals." In
%   J. Blauert, editor, The technology of binaural listening, chapter 7.
%   Springer, Berlin, Heidelberg, New York NY, 2013.
%   url: http://link.springer.com/chapter/10.1007%2F978-3-642-37762-4_7
%
%   [2] E. Georganti, T. May, S. van de Par, and J. Mourjopoulos, "Sound
%   Source Distance Estimation in Rooms based on Statistical Properties of
%   Binaural Signals", IEEE Transactions on Audio, Speech and Language 
%   Processing, Vol. 21 (8), Aug. 2013.
%   url: http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6508870
%
%   ---------
%
%   REQUIREMENTS: 
%
%   1) Download audio files from 
%      http://sourceforge.net/projects/audiofordistanceestimation/files/wavFilesGeorganti.zip/download
%      (copy/paste the above link at your browser)
%
%   2) Unzip folder (wavFilesGeorganti.zip)
%
%   3) Move folder in amtoolbox/signals/
%
%
%   ---------
%
%   EXAMPLE:
%
%   To display Fig.9 [1] use :
%
%     exp_georganti2013('fig9');
%
%   To display Fig.10 [1] use :
%
%     exp_georganti2013('fig10');
%
%
%   Url: http://amtoolbox.sourceforge.net/doc/experiments/exp_georganti2013.php

%   Copyright (C) 2009-2014 Peter L. S??ndergaard and Piotr Majdak.
%   This file is part of AMToolbox version 0.9.5
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%
%   Developed with Matlab 7.1.0.584 (R2010b)      v.0.1
%   Updated   with Matlab R2011b    (7.13.0.564)  v.1.0
%
%   Please send bug reports to:
%
%   Author  :  Eleftheria Georganti
%              Postdoctoral Researcher
%              Experimental Audiology, ENT
%              University Hospital of Zurich/University of Zurich
%              Zurich, Switzerland
%              eleftheria.georganti@uzh.ch

%   History :
%   v.0.1   2013/01/21
%   v.1.0   2014/06/10 Link to wavfiles & documentation
% 
%
%*************************************************************************



% SELECT SMALL OR LARGE ROOM BY COMMENTING/UNCOMMENTING
definput.flags.type = {'missingflag','fig9','fig10'};

% Parse input options
[flags]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

if flags.do_fig9
    num = 1;      % Small Room -  RT = 0.15 sec
    prefix = 'smallRoom';
elseif flags.do_fig10
    num = 2;      % Large Room -  RT = 0.9 sec
    prefix = 'largeRoom';
end

% Sampling frequency
P.fs = 44100;

% Frame size in seconds
P.timeFr = 1;

% Frame size in samples
P.sampleFr = P.fs * P.timeFr;

% Overlap
P.hop = P.sampleFr/2;

% FFT points
P.nFFT = P.sampleFr;

% Frequency index
P.freq = (P.fs/P.nFFT)*(0:(P.nFFT/2-1));

% Define frequency range for the BSDM STD feature calculation
P.fmin = 20;    % lower frequency in Hz - default value
P.fmax = 2300; % upper frequency in Hz - default value

nDist=dir(fullfile(amtbasepath,'signals',['exp_georganti2013_' prefix '*.mat']));

fmin_id = min(find((P.freq>P.fmin)));
fmax_id = min(find((P.freq>P.fmax)));


idx = 1;

for ww = 1:length(nDist)
    
    pathWav = nDist(ww).name;
    signalPre = load(pathWav);
    signal = signalPre.signal;
    
    for kk = 1:P.hop:length(signal)-P.hop
        
        
        % Calculate magnitude spectrums in dB of the left & right signals
        leftFFT  = 20*log10(abs(fft(signal(kk:kk+P.hop-1,1))));
        rightFFT = 20*log10(abs(fft(signal(kk:kk+P.hop-1,2))));
        
        % Subtract the magnitude spectrums
        specDIF  = leftFFT(1:end/2)-rightFFT(1:end/2);
        
        % Calculate the differential standard deviation for the
        % frequency range of interest
        difSTD(ww,idx) = std(specDIF(fmin_id:fmax_id));
        
        clear leftFFT rightFFT specDIFF
        
        idx = idx+1;
        
    end
    
    idx = 1;
end



% Plot figures

if num == 1
    figure;
    subplot(1,3,[1 2])
    plot(difSTD(1,10:210),'k','LineWidth',2)
    hold on
    plot(difSTD(2,10:210),'Color',[0.6 0.6 0.6],'LineWidth',2)
    plot(difSTD(3,10:210),'--k','LineWidth',2)
    xlabel('Time (sec)','FontSize',12,'FontWeight','bold')
    ylabel('BSMD STD','FontSize',12,'FontWeight','bold')
    title ('Small Room -  RT = 0.15 sec','FontSize',16,'FontWeight','bold')
    ylim([3 8])
    xlim([0 200])
    legend('0.5m','1m','1.5m','Location','NorthWest')
    set(gca,'FontSize',10)
    grid on
    box on
    
    subplot(1,3,3)
    [a1,b1] = hist(difSTD(1,10:210),[3:0.1:8]);
    [a2,b2] = hist(difSTD(2,10:210),[3:0.1:8]);
    [a3,b3] = hist(difSTD(3,10:210),[3:0.1:8]);
    plot(a1/201,b1,'LineWidth',2,'Color','k'); hold on
    set(gca,'FontSize',10)
    plot(a2/201,b2,'Color',[0.6 0.6 0.6],'LineWidth',2)
    plot(a3/201,b3,'LineWidth',2,'Color','k','LineStyle','--')
    xlabel('Frequency of occurence (%)','FontSize',10,'FontWeight','bold')
    legend('0.5m','1m','1.5m','Location','NorthEast')
    box on
    xlim([0 0.4])
    ylim([3 8])
    box on
    set(gca,'YTick',[])
    grid on
    
%     print -depsc2 smallSpeechBSMD_STD  %Creates .eps file with the figure
    
    
elseif num == 2 
    
    figure;
    subplot(1,3,[1 2])
    plot(difSTD(1,10:210),'k','LineWidth',2)
    hold on
    plot(difSTD(2,10:210),'Color',[0.6 0.6 0.6],'LineWidth',2)
    plot(difSTD(3,10:210),'--k','LineWidth',2)
    xlabel('Time (sec)','FontSize',12,'FontWeight','bold')
    ylabel('BSMD STD','FontSize',12,'FontWeight','bold')
    title ('Large Room -  RT = 0.9 sec','FontSize',16,'FontWeight','bold')
    ylim([2 8.5])
    xlim([0 200])
    legend('1m','2m','3m','Location','NorthWest')
    set(gca,'FontSize',10)
    grid on
    box on
    
    subplot(1,3,3)
    [a1,b1] = hist(difSTD(1,10:210),[3:0.1:8]);
    [a2,b2] = hist(difSTD(2,10:210),[3:0.1:8]);
    [a3,b3] = hist(difSTD(3,10:210),[3:0.1:8]);
    plot(a1/201,b1,'LineWidth',2,'Color','k'); hold on
    set(gca,'FontSize',10)
    plot(a2/201,b2,'Color',[0.6 0.6 0.6],'LineWidth',2)
    plot(a3/201,b3,'LineWidth',2,'Color','k','LineStyle','--')
    xlabel('Frequency of occurence (%)','FontSize',10,'FontWeight','bold')
    legend('1m','2m','3m','Location','NorthEast')
    box on
    xlim([0 0.4])
    ylim([2 8.5])
    box on
    set(gca,'YTick',[])
    grid on
    
end




