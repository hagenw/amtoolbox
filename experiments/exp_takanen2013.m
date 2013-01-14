function output = exp_takanen2013(varargin)
%EXP_TAKANEN2013 Figures from Takanen, Santala, Pulkki (2013a,2013b)
%Usage: output = exp_takanen2013(flag)
%
%   EXP_TAKANEN2013(flag) reproduces the figure given by flag either from
%   the Takanen et al. (2013) book chapter or the Takanen et al. (2013)
%   manuscript. The format of its output depends on the chosen figure.
%   Optionally, pre-computed cochlear model outputs for the different
%   scenarios can be applied to significantly reduce the required
%   computation time. The pre-computed cochlear model outputs can be 
%   obtained from the authors.
%  
%   The following flags can be specified:
%
%     'binsig'      : use binaural input signals in the computation. This 
%                     is the default.
%
%     'cochlea'     : use pre-computed cochlea model outputs in the
%                     computation to reduce computation time. 
%
%     'fig7book'    : Figure 7 from the book chapter. Binaural activity 
%                     maps obtained with the model for an off-sweet-spot 
%                     listening scenario with different audio coding 
%                     techniques.
%
%     'fig8book'    : Figure 8 from the book chapter. Activation
%                     distributions obtained with the model for (a) the 
%                     reference scenario of incoherent pink noise emitted 
%                     from twelve azimuth directions, and (b)-(d) the 
%                     reproduction of such a scenario with an eight-channel
%                     loudspeaker system employing signals obtained with
%                     different audio coding techniques. Additionally, the
%                     the distributions when DirAC is used in audio coding
%                     of 5.0 surround signal having incoherent pink noise
%                     in each channel with (e) the straightforward method 
%                     and (f) the even-layout method.
%
%     'fig6art'     : Figure 6 from the manuscript. Binaural activity maps 
%                     for four binaural listening scenarios, namely (a)
%                     HRTF-processed pink noise, (b) pink noise with ITD, 
%                     (c) anti-phasic sinusoidal sweep, and (d) band-
%                     limited noise centered around 500 Hz with an ITD of
%                     1.5 ms.
%
%     'fig7art'     : Figure 7 from the manuscript. Binaural activity maps 
%                     for four binaural listening scenarios, namely (a) 
%                     $S_\pi N_0$ with different signal-to-noise ratios, 
%                     (b) binaural interference, (c) precedence effect, and
%                     (d) binaural room impulse response.
%
%   If no flag is given, the function will print the list of valid flags.
%
%   Examples:
%   ---------
%
%   To display Figure 7 from the book chapter using pre-computed cochlea
%   model outputs use:
%
%     exp_takanen2013('fig7book','cochlea');
%
%   To display Figure 8 from the book chapter using pre-computed cochlea
%   model outputs use:
%
%     exp_takanen2013('fig8book','cochlea');
%
%   To display Figure 6 from the manuscript using pre-computed cochlea 
%   model outputs use:
%
%     exp_takanen2013('fig6art','cochlea');
%
%   To display Figure 6 from the manuscript using pre-computed cochlea 
%   model outputs use:
%
%     exp_takanen2013('fig7art','cochlea');
%
%   References: takanen2013a takanen2013b
%
%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland

definput.import={'amtredofile'};
definput.flags.type={'missingflag','fig7book','fig8book','fig6art','fig7art'};

definput.flags.dataType={'binsig','cochlea'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% Setting of parameters
fs = 48000;
printFigs = 0;
printMap =0;
compType =1;
h = figure;
%% Figure 7 from the book chapter
if flags.do_fig7book
    % if the user wishes to compute the cochlear model outputs, binaural
    % input signals are used
    if flags.do_binsig
        load exp_takanen2013fig7bookbinsignals.mat -mat
        for ind=1:length(tests)
            % compute the binaural activity map with the model
            output = takanen2013(tests(ind).insig,fs,compType,printFigs,printMap);
            nXBins= length(output.levels)*(size(output.colorMtrx,1)-1);
            dim = size(output.activityMap);
            output.colorGains(output.colorGains>1) =1;
            outputMtrx = zeros(dim(1),nXBins,3);
            for colorInd=1:size(output.colorMtrx,1)
                temp = find((output.activityMap==(colorInd-1))==1);
                outputMtrx(temp) = output.colorGains(temp)*output.colorMtrx(colorInd,1);
                outputMtrx(temp+dim(1)*nXBins) = output.colorGains(temp)*output.colorMtrx(colorInd,2);
                outputMtrx(temp+2*dim(1)*nXBins) = output.colorGains(temp)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(3,2,ind);imagesc(output.levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
            title(tests(ind).case);
            set(gca,'YTick',.0:.5:2.5);
            set(gca,'YTickLabel',2.5:-0.5:0);
            set(gca,'Xtick',-1:0.4:1);
            xlabel('Activation location');
            ylabel('Time [s]');
        end
        
    end
    %otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea
        load exp_takanen2013fig7bookcochleadata -mat
        for ind=1:length(tests)
            % compute the binaural activity map with the model
            output = takanen2013(tests(ind).cochlear,fs,compType,printFigs,printMap);
            nXBins= length(output.levels)*(size(output.colorMtrx,1)-1);
            dim = size(output.activityMap);
            output.colorGains(output.colorGains>1) =1;
            outputMtrx = zeros(dim(1),nXBins,3);
            for colorInd=1:size(output.colorMtrx,1)
                temp = find((output.activityMap==(colorInd-1))==1);
                outputMtrx(temp) = output.colorGains(temp)*output.colorMtrx(colorInd,1);
                outputMtrx(temp+dim(1)*nXBins) = output.colorGains(temp)*output.colorMtrx(colorInd,2);
                outputMtrx(temp+2*dim(1)*nXBins) = output.colorGains(temp)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(3,2,ind);imagesc(output.levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
            title(tests(ind).case);
            set(gca,'YTick',.0:.5:2.5);
            set(gca,'YTickLabel',2.5:-0.5:0);
            set(gca,'Xtick',-1:0.4:1);
            ylabel('Time [s]');
            xlabel('Activation location');
        end
    end
end
%% Figure 8 from the book chapter
if flags.do_fig8book
    probDist = zeros(6,19);
    % if the user wishes to compute the cochlear model outputs, binaural
    % input signals are used
    if flags.do_binsig
        load exp_takanen2013fig8bookbinsignals.mat -mat
        for ind=1:length(tests)
            % compute the binaural activity map with the model
            output = takanen2013(tests(ind).insig,fs,compType,printFigs,printMap);
            for i=1:6
                probDist(i,:) = sum(output.colorGains(:,i:6:end));
            end
            temp = probDist./(max(probDist,[],2)*ones(1,size(probDist,2)));
            outputMtrx = zeros(size(temp,1),size(temp,2),3);
            for colorInd=2:size(output.colorMtrx,1)
                outputMtrx(colorInd-1,:,1) = temp(colorInd-1,:)*output.colorMtrx(colorInd,1);
                outputMtrx(colorInd-1,:,2) = temp(colorInd-1,:)*output.colorMtrx(colorInd,2);
                outputMtrx(colorInd-1,:,3) = temp(colorInd-1,:)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(3,2,ind);imagesc(output.levels./90,6:-1:1,outputMtrx);
            title(tests(ind).case);
            set(gca,'YTick',1:6);
            set(gca,'YTickLabel',6:-1:1);
            set(gca,'Xtick',-1:0.4:1);
            ylabel('Frequency area');
            xlabel('Distribution of activation');
        end
    end
    %otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea
        load exp_takanen2013fig8bookcochleadata -mat
        for ind=1:length(tests)
            % compute the binaural activity map with the model
            output = takanen2013(tests(ind).cochlear,fs,compType,printFigs,printMap);
            for i=1:6
                probDist(i,:) = sum(output.colorGains(:,i:6:end));
            end
            temp = probDist./(max(probDist,[],2)*ones(1,size(probDist,2)));
            outputMtrx = zeros(size(temp,1),size(temp,2),3);
            for colorInd=2:size(output.colorMtrx,1)
                outputMtrx(colorInd-1,:,1) = temp(colorInd-1,:)*output.colorMtrx(colorInd,1);
                outputMtrx(colorInd-1,:,2) = temp(colorInd-1,:)*output.colorMtrx(colorInd,2);
                outputMtrx(colorInd-1,:,3) = temp(colorInd-1,:)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(3,2,ind);imagesc(output.levels./90,6:-1:1,outputMtrx);
            title(tests(ind).case);
            set(gca,'YTick',1:6);
            set(gca,'YTickLabel',6:-1:1);
            set(gca,'Xtick',-1:0.4:1);
            ylabel('Frequency area');
            xlabel('Distribution of activation');
        end
    end
end
%% Figure 6 from the article
if flags.do_fig6art
    % if the user wishes to compute the cochlear model outputs, binaural
    % input signals are used
    if flags.do_binsig
        load exp_takanen2013fig6artbinsignals.mat -mat
        for ind=1:length(tests)
            activityMap = [];
            gains = [];
            %some scenarios consist of multiple test cases that are
            %processed separately
            for caseInd=1:length(tests(ind).binSignals)
                % compute the binaural activity map with the model
                output = takanen2013(tests(ind).binSignals(caseInd).insig,fs,compType,printFigs,printMap);
                %concatenate the separate activity maps into one map
                activityMap = [activityMap;output.activityMap];
                gains = [gains;output.colorGains];
            end
            %the anti-phasic sweep contains also frequencies below the
            %frequency range of the model. Hence, the first 0.5 s of the
            %activity map are removed
            if(strcmp('Anti-phasic sinusoidal sweep',tests(ind).scenario)==1)
                activityMap = activityMap(0.5*fs+1:end,:);
                gains = gains(0.5*fs+1:end,:);
            end
            nXBins= length(output.levels)*(size(output.colorMtrx,1)-1);
            dim = size(activityMap);
            gains(gains>1) =1;
            outputMtrx = zeros(dim(1),nXBins,3);
            for colorInd=1:size(output.colorMtrx,1)
                temp = find((activityMap==(colorInd-1))==1);
                outputMtrx(temp) = gains(temp)*output.colorMtrx(colorInd,1);
                outputMtrx(temp+dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,2);
                outputMtrx(temp+2*dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(2,2,ind);imagesc(output.levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
            title(tests(ind).scenario);
            set(gca,'YTick',tests(ind).ytickPos);
            set(gca,'YTickLabel',tests(ind).ytickLab(end:-1:1));
            set(gca,'Xtick',-1:0.4:1);
            ylabel(tests(ind).ylab);
            xlabel('Activation location');
        end
    end
    %otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea
        load exp_takanen2013fig6artcochleadata -mat
        for ind=1:length(tests)
            activityMap = [];
            gains = [];
            %some scenarios consist of multiple test cases that are
            %processed separately
            for caseInd=1:length(tests(ind).cochlearData)
                % compute the binaural activity map with the model
                output = takanen2013(tests(ind).cochlearData(caseInd).cochlear,fs,compType,printFigs,printMap);
                %concatenate the separate activity maps into one map
                activityMap = [activityMap;output.activityMap];
                gains = [gains;output.colorGains];
            end
            %the anti-phasic sweep contains also frequencies below the
            %frequency range of the model. Hence, the first 0.5 s of the
            %activity map are removed
            if(strcmp('Anti-phasic sinusoidal sweep',tests(ind).scenario)==1)
                activityMap = activityMap(0.5*fs+1:end,:);
                gains = gains(0.5*fs+1:end,:);
            end
            nXBins= length(output.levels)*(size(output.colorMtrx,1)-1);
            dim = size(activityMap);
            gains(gains>1) =1;
            outputMtrx = zeros(dim(1),nXBins,3);
            for colorInd=1:size(output.colorMtrx,1)
                temp = find((activityMap==(colorInd-1))==1);
                outputMtrx(temp) = gains(temp)*output.colorMtrx(colorInd,1);
                outputMtrx(temp+dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,2);
                outputMtrx(temp+2*dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(2,2,ind);imagesc(output.levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
            title(tests(ind).scenario);
            set(gca,'YTick',tests(ind).ytickPos);
            set(gca,'YTickLabel',tests(ind).ytickLab(end:-1:1));
            set(gca,'Xtick',-1:0.4:1);
            ylabel(tests(ind).ylab);
            xlabel('Activation location');
        end
    end
end
%% Figure 7 from the article
if flags.do_fig7art
    % if the user wishes to compute the cochlear model outputs, binaural
    % input signals are used
    if flags.do_binsig
        load exp_takanen2013fig7artbinsignals.mat -mat
        for ind=1:length(tests)
            activityMap = [];
            gains = [];
            %some scenarios consist of multiple test cases that are
            %processed separately
            for caseInd=1:length(tests(ind).binSignals)
                % compute the binaural activity map with the model
                output = takanen2013(tests(ind).binSignals(caseInd).insig,fs,compType,printFigs,printMap);
                %concatenate the separate activity maps into one map
                activityMap = [activityMap;output.activityMap];
                gains = [gains;output.colorGains];
            end
            nXBins= length(output.levels)*(size(output.colorMtrx,1)-1);
            %in order to better visualize the clicks in the precedence
            %effect scenario, most of the silent parts of the signal
            %are removed
            if(strcmp('Precedence effect',tests(ind).scenario)==1)
                activityMap = activityMap([1500:3700 4500:6700 7500:9700 10500:12700 13500:15700 16500:18700 20200:22400],:);
                gains = gains([1500:3700 4500:6700 7500:9700 10500:12700 13500:15700 16500:18700 20200:22400],:);
                gains = 2*gains;
            end
            dim = size(activityMap);
            gains(gains>1) =1;
            outputMtrx = zeros(dim(1),nXBins,3);
            for colorInd=1:size(output.colorMtrx,1)
                temp = find((activityMap==(colorInd-1))==1);
                outputMtrx(temp) = gains(temp)*output.colorMtrx(colorInd,1);
                outputMtrx(temp+dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,2);
                outputMtrx(temp+2*dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(2,2,ind);imagesc(output.levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
            title(tests(ind).scenario);
            set(gca,'YTick',tests(ind).ytickPos);
            set(gca,'YTickLabel',tests(ind).ytickLab);
            set(gca,'Xtick',-1:0.4:1);
            ylabel(tests(ind).ylab);
            xlabel('Activation location');
        end
    end
    %otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea
        load exp_takanen2013fig7artcochleadata -mat
        for ind=1:length(tests)
            activityMap = [];
            gains = [];
            %some scenarios consist of multiple test cases that are
            %processed separately
            for caseInd=1:length(tests(ind).cochlearData)
                % compute the binaural activity map with the model
                output = takanen2013(tests(ind).cochlearData(caseInd).cochlear,fs,compType,printFigs,printMap);
                %concatenate the separate activity maps into one map
                activityMap = [activityMap;output.activityMap];
                gains = [gains;output.colorGains];
            end
            nXBins= length(output.levels)*(size(output.colorMtrx,1)-1);
            %in order to better visualize the clicks in the precedence
            %effect scenario, most of the silent parts of the signal
            %are removed
            if(strcmp('Precedence effect',tests(ind).scenario)==1)
                activityMap = activityMap([1500:3700 4500:6700 7500:9700 10500:12700 13500:15700 16500:18700 20200:22400],:);
                gains = gains([1500:3700 4500:6700 7500:9700 10500:12700 13500:15700 16500:18700 20200:22400],:);
                gains = 2*gains;
            end
            dim = size(activityMap);
            gains(gains>1) =1;
            outputMtrx = zeros(dim(1),nXBins,3);
            for colorInd=1:size(output.colorMtrx,1)
                temp = find((activityMap==(colorInd-1))==1);
                outputMtrx(temp) = gains(temp)*output.colorMtrx(colorInd,1);
                outputMtrx(temp+dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,2);
                outputMtrx(temp+2*dim(1)*nXBins) = gains(temp)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(2,2,ind);imagesc(output.levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
            title(tests(ind).scenario);
            %suplementary data is plotted in the binaural
            %interference scenario
            if(isempty(tests(ind).cases)==0)
                text(-.8,2,tests(ind).cases(1),'horizontalAlignment','left','verticalAlignment','top','Rotation',0,'color','w');
                text(-.8,.2,tests(ind).cases(2),'horizontalAlignment','left','verticalAlignment','top','Rotation',0,'color','w');
                hold on;
                plot([-1 1],[1.9 1.9],'LineStyle','--','Color','w')
                hold off;
            end
            set(gca,'YTick',tests(ind).ytickPos);
            set(gca,'YTickLabel',tests(ind).ytickLab);
            ylabel(tests(ind).ylab);
            set(gca,'Xtick',-1:0.4:1);
            xlabel('Activation location');
        end
    end
end
output = g;