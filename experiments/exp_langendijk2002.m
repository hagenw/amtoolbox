function output = exp_langendijk2002(varargin)
%EXP_LANGENDIJK2002  Demo of the localization model of Langendijk & Bronkhorst 2002
%
%   Validation of Langendijk et al. (2002)
%
%   This script generates figures showing the results of the localization 
%   model based on the simulation shown in the related paper of Langendijk &
%   Bronkhorst 2002. 
%
%   The following flags can be specified;
%
%-    'plot' - plot the output of the experiment. This is the default.
%
%-    'noplot' - don't plot, only return data.
%
%-    'fig7' - listener P6
%
%-    'fig9' - listener P3
%
%   You can choose between two of his listeners P3 and P6. The required
%   data (DTF data and response patterns) will be provided by precalculated
%   mat-files due to high computing time (optionally data can be calculated
%   by using the data_langendijk2002 function). 
%
%
%   subfigure 1 Baseline condition
%   subfigure 2 2-octave condition (4-16kHz)
%   subfigure 3 1-octave condition (low 4-8kHz)
%   subfigure 4 1-octave condition (middle 5.7-11.3kHz)
%   subfigure 5 1-octave condition (high 8-16kHz)
%
%     Above-named subfigures are showing the probability density function (pdf) 
%     and actual responses(°) for chosen listener as a function of target 
%     position for different conditions. The shading of each cell codes the 
%     probability density (light/dark is high/low probability)
%
%   subfigure 6 Likelihood statistics
%
%     subfigure 6 shows the likelihood statistics for the actual responses
%     (bars), the average and the 99% confidence interval of the expected
%     likelihood (dots and bars). See the paper for further details.
%
%   The output are the pdfs for the baseline condition.
%
%   See also: langendijk, likelilangendijk, plotlangendijk, plotlikelilangendijk

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ------ Check input options --------------------------------------------

  definput.flags.type = {'missingflag','fig7','fig9'};
  definput.flags.plot = {'plot','noplot'};

  % Parse input options
  [flags,keyvals]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig7
    listener='P6';  % ID of listener (P3 or P6)
elseif flags.do_fig9
    listener='P3';
end

fs = 48000;     % sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['langendijk2002-' listener]); 
% loads hM data for all conditions which were calculated as follows:
% temp=data_langendijk2002([listener '-dtf']);
% pol=temp(1,:);
% med=temp(2:end,:);
% temp=data_langendijk2002([listener '-b']);
% targetb=temp(1,:); responseb=temp(2,:);
% medir=gr2ir(med,'b',fs);
% temp=data_langendijk2002([listener '-2o']);
% target2o=temp(1,:); response2o=temp(2,:);
% medir2o=gr2ir(med,'2o',fs);
% temp=data_langendijk2002([listener '-1ol']);
% target1ol=temp(1,:); response1ol=temp(2,:);
% medir1ol=gr2ir(med,'1ol',fs);
% temp=data_langendijk2002([listener '-1om']);
% target1om=temp(1,:); response1om=temp(2,:);
% medir1om=gr2ir(med,'1om',fs);
% temp=data_langendijk2002([listener '-1oh']);
% target1oh=temp(1,:); response1oh=temp(2,:);
% medir1oh=gr2ir(med,'1oh',fs);

% pdf calcualtion
h = waitbar(0,'Please wait...');
pb  = langendijk( medir   ,medir); % baseline
waitbar(1/5)
p2o = langendijk( medir2o ,medir); % 2-oct (4-16kHz)
waitbar(2/5)
p1ol= langendijk( medir1ol,medir); % 1-oct (low:4-8kHz)
waitbar(3/5)
p1om= langendijk( medir1om,medir); % 1-oct (middle:5.7-11.3kHz)
waitbar(4/5)
p1oh= langendijk( medir1oh,medir); % 1-oct (high:8-16kHz)
waitbar(5/5)

% likelihood estimations
la=zeros(5,1);le=zeros(5,1);ci=zeros(5,2);
idb=1:2:length(targetb); % in order to get comparable likelihoods
[la(1),le(1),ci(1,:)] = likelilangendijk( pb,pol,pol,targetb(idb),responseb(idb) );
[la(2),le(2),ci(2,:)] = likelilangendijk( p2o,pol,pol,targetc,response2o );
[la(3),le(3),ci(3,:)] = likelilangendijk( p1ol,pol,pol,targetc,response1ol );
[la(4),le(4),ci(4,:)] = likelilangendijk( p1om,pol,pol,targetc,response1om );
[la(5),le(5),ci(5,:)] = likelilangendijk( p1oh,pol,pol,targetc,response1oh );
close(h)

output = pb;

if flags.do_plot
    figure
    clf
    
    % pdf plots with actual responses
    subplot(2,3,1)
    plotlangendijk(pb,pol,pol,[listener '; ' 'baseline']);
    hold on; h=plot( targetb, responseb, 'ko'); set(h,'MarkerFaceColor','w')
    subplot(2,3,2)
    plotlangendijk(p2o,pol,pol,[listener '; ' '2-oct (4-16kHz)']);
    hold on; h=plot( targetc, response2o, 'ko'); set(h,'MarkerFaceColor','w')
    subplot(2,3,3)
    plotlangendijk(p1ol,pol,pol,[listener '; ' '1-oct (low: 4-8kHz)']);
    hold on; h=plot( targetc, response1ol, 'ko'); set(h,'MarkerFaceColor','w')
    subplot(2,3,4)
    plotlangendijk(p1om,pol,pol,[listener '; ' '1-oct (middle: 5.7-11.3kHz)']);
    hold on; h=plot( targetc, response1om, 'ko'); set(h,'MarkerFaceColor','w')
    subplot(2,3,5)
    plotlangendijk(p1oh,pol,pol,[listener '; ' '1-oct (high: 8-16kHz)']);
    hold on; h=plot( targetc, response1oh, 'ko'); set(h,'MarkerFaceColor','w')
    
    % likelihood statistic
    subplot(2,3,6)
    plotlikelilangendijk(la,le,ci)
end


function [medir]=gr2ir(med,cond,fs)
% GR2IR converts given gain responses MED (in dB) to impulse responses
% MEDIR; furthermore several conditions according to langendijk et al.
% (2002) can be defined
% Usage:            [medir]=gr2ir(med,cond,fs)
% Input arguments:
%       med:        gain responses (in dB)
%       cond:       condition, 
%                   possibilities:  baseline    'b'
%                                   2 octaves   '2o'
%                                   1 oct (low) '1ol'
%                                   1 oct (mid) '1om'
%                                   1 oct (high)'1oh'
%       fs:      	sampling frequency
% Output arguments:
%       medir:   	impulse responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default settings
if ~exist('cond','var')
    cond='b';
end
if ~exist('fs','var')
    fs=48000;
end

imp=zeros(240,1);imp(1)=1;
n=256;
len=240;
frq=logspace(log10(2000),log10(16000),size(med,1));

% frequency indices
f1=find(frq>=4000,1);
f2=find(frq>=5700,1);
f3=find(frq>=8000,1);
f4=find(frq>=11300,1);

switch cond
    case 'b' % baseline
        medir=zeros(n,size(med,2));
        for ii=1:size(med,2)
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '2o' % 2 oct (4-16kHz)
        medir=zeros(n,size(med,2));
        med2o=med;
        for ii=1:size(med,2)
            med2o(f1:end,ii)=mean(med(f1:end,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med2o(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1ol' % 1 oct (low:4-8kHz)
        medir=zeros(n,size(med,2));
        med1ol=med;
        for ii=1:size(med,2)
            med1ol(f1:f3,ii)=mean(med(f1:f3,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1ol(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1om' % 1 oct (middle:5.7-11.3kHz)
        medir=zeros(n,size(med,2));
        med1om=med;
        for ii=1:size(med,2)
            med1om(f2:f4,ii)=mean(med(f2:f4,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1om(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1oh' % 1 oct (high:8-16kHz)
        medir=zeros(n,size(med,2));
        med1oh=med;
        for ii=1:size(med,2)
            med1oh(f3:end,ii)=mean(med(f3:end,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1oh(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end
end
end
