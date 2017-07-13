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
%   The figures show the computed distance-dependent feature
%   BSMD-STD (Binaural spectral-magnitude difference standard deviation)
%   for audio signals measured at different/receiver distances.
%
%   The BSMD-STD may be derived from any dual-channel signal
%   (binaural/stereo recordings).
%
%   The BSMD-STD is related to the standard deviation of the
%   magnitude spectrum of room impulse response, which is known to depend on 
%   the source/receiver distance. See Georganti et al. (2013) for more
%   information on applications.
%
%
%   Examples:
%   ---------
%
%   To display Fig. 9 use :::
%
%     exp_georganti2013('fig9');
%
%   To display Fig. 10 use :::
%
%     exp_georganti2013('fig10');
%
%   See also: georganti2013
%
%   References: georganti2013 georganti2013application
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
%   v.1.1   2015/01/20 adapted to AMT (Piotr Majdak)
% 
%
%*************************************************************************

definput.flags.type = {'missingflag','fig9','fig10'};

% Parse input options
flags  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

  % Set model parameters
P.fs = 44100; % Sampling frequency
P.timeFr = 1; % Frame size in seconds
P.fmin = 20;    % lower frequency in Hz for the BSDM STD calculation
P.fmax = 2300; % upper frequency in Hz for the BSDM STD calculation


if flags.do_fig9 % Small Room -  RT = 0.15 sec
      % load signals        
    nDist(1).name='smallRoomSpeech_0_5m.mat';
    nDist(2).name='smallRoomSpeech_1_0m.mat';
    nDist(3).name='smallRoomSpeech_1_5m.mat';
      % calculate 
    for ww = 1:length(nDist)   
        signalPre = amt_load('georganti2013',nDist(ww).name);    
        BSMDSTD(ww,:) = georganti2013(signalPre.signal,P);
    end
      % plot BSMD-STDs
    figure;
    subplot(1,3,[1 2])
    plot(BSMDSTD(1,10:210),'k','LineWidth',2)
    hold on
    plot(BSMDSTD(2,10:210),'Color',[0.6 0.6 0.6],'LineWidth',2)
    plot(BSMDSTD(3,10:210),'--k','LineWidth',2)
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
    [a1,b1] = hist(BSMDSTD(1,10:210),3:0.1:8);
    [a2,b2] = hist(BSMDSTD(2,10:210),3:0.1:8);
    [a3,b3] = hist(BSMDSTD(3,10:210),3:0.1:8);
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
    
end

if flags.do_fig10 % Large Room -  RT = 0.9 sec
      % load signals        
    nDist(1).name='largeRoomSpeech_1m.mat';
    nDist(2).name='largeRoomSpeech_2m.mat';
    nDist(3).name='largeRoomSpeech_3m.mat';
      % calculate 
    for ww = 1:length(nDist)   
        signalPre = amt_load('georganti2013',nDist(ww).name);    
        BSMDSTD(ww,:) = georganti2013(signalPre.signal,P);
    end
      % plot BSMD-STDs
    figure;
    subplot(1,3,[1 2])
    plot(BSMDSTD(1,10:210),'k','LineWidth',2)
    hold on
    plot(BSMDSTD(2,10:210),'Color',[0.6 0.6 0.6],'LineWidth',2)
    plot(BSMDSTD(3,10:210),'--k','LineWidth',2)
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
    [a1,b1] = hist(BSMDSTD(1,10:210),[3:0.1:8]);
    [a2,b2] = hist(BSMDSTD(2,10:210),[3:0.1:8]);
    [a3,b3] = hist(BSMDSTD(3,10:210),[3:0.1:8]);
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

