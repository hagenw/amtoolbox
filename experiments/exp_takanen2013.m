function output = exp_takanen2013(varargin)
%EXP_TAKANEN2013 Figures from Takanen, Santala, Pulkki (2013a,2013b)
%   Usage: output = exp_takanen2013(flag)
%
%   `exp_takanen2013(flag)` reproduces the figure given by *flag* either from
%   the Takanen et al. (2013) book chapter or the Takanen et al. (2014)
%   manuscript. The format of its output depends on the chosen figure.
%   Optionally, pre-computed cochlear model outputs for the different
%   scenarios can be applied to significantly reduce the required
%   computation time. The pre-computed cochlear model outputs can be 
%   obtained from the authors.
%  
%   The following flags can be specified:
%
%     'binsig'       use binaural input signals in the computation. This 
%                    is the default.
%
%     'cochlea'      use pre-computed cochlea model outputs in the
%                    computation to reduce computation time. 
%
%     'fig8'     Figure 8 from the book chapter Takanen et al. (2013). Binaural activity 
%                    maps obtained with the model for an off-sweet-spot 
%                    listening scenario with different audio coding 
%                    techniques.
%
%     'fig9'     Figure 9 from the book chapter Takanen et al. (2013). Activation
%                    distributions obtained with the model for (a) the 
%                    reference scenario of incoherent pink noise emitted 
%                    from twelve azimuth directions, and (b)-(d) the 
%                    reproduction of such a scenario with an eight-channel
%                    loudspeaker system employing signals obtained with
%                    different audio coding techniques. Additionally, the
%                    the distributions when DirAC is used in audio coding
%                    of 5.0 surround signal having incoherent pink noise
%                    in each channel with (e) the straightforward method 
%                    and (f) the even-layout method.
%
%     'fig7_takanen2014' Figure 7 from the article Takanen et al. (2014). Binaural activity maps 
%                    for four binaural listening scenarios, namely (a)
%                    HRTF-processed pink noise, (b) pink noise with ITD, 
%                    (c) anti-phasic sinusoidal sweep, and (d) band-
%                    limited noise centered around 500 Hz with an ITD of
%                    1.5 ms.
%
%     'fig8_takanen2014' Figure 8 from the article Takanen et al. (2014). Binaural activity maps 
%                    for four binaural listening scenarios, namely (a) 
%                    $S_\pi N_0$ with different signal-to-noise ratios, 
%                    (b) binaural interference, (c) precedence effect, and
%                    (d) binaural room impulse response.
%
%   If no flag is given, the function will print the list of valid flags.
%
%   Requirements and installation: 
%   1) Functioning model verhulst2012 (see the corresponding requirements)
%
%   2) Data from www.acoustics.hut.fi/publications/papers/AMTool2013-bam/ in amtbase/signals
%
%   3) at least 3 GB of RAM
%
%   Examples:
%   ---------
%
%   To display Figure 8 from the book chapter Takanen et al. (2013) using pre-computed cochlea
%   model outputs use:::
%
%     exp_takanen2013('fig8','cochlea');
%
%   To display Figure 9 from the book chapter Takanen et al. (2013) using pre-computed cochlea
%   model outputs use:::
%
%     exp_takanen2013('fig9','cochlea');
%
%   To display Figure 7 from the article Takanen et al. (2014) using pre-computed cochlea 
%   model outputs use:::
%
%     exp_takanen2013('fig7_takanen2014','cochlea');
%
%   To display Figure 8 the article Takanen et al. (2014) using pre-computed cochlea 
%   model outputs use:::
%
%     exp_takanen2013('fig8_takanen2014','cochlea');
%
%   References: takanen2013 takanen2014

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland

definput.import={'amtredofile'};
definput.flags.type={'missingflag','fig8','fig9','fig7_takanen2014','fig8_takanen2014'};

definput.flags.dataType={'cochlea','binsig'};

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
%% Figure 8 from the book chapter
if flags.do_fig8
    % if the user wishes to compute the cochlear model outputs, binaural
    % input signals are used
    if flags.do_binsig
        data=safe_load('exp_takanen2013fig8bookbinsignals.mat');
        for ind=1:length(data.tests)
            data=safe_load('exp_takanen2013fig8bookbinsignals.mat');
            insig=data.tests(ind).insig;
            tit=data.tests(ind).case;
            clear data
            % compute the binaural activity map with the model
            output = takanen2013(insig,fs,compType,printFigs,printMap);
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
            clear output outputMtrx
            title(tit);
            set(gca,'YTick',.0:.5:2.5);
            set(gca,'YTickLabel',2.5:-0.5:0);
            set(gca,'Xtick',-1:0.4:1);
            xlabel('Activation location');
            ylabel('Time [s]');
        end
        
    end
    %otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea
        data=safe_load('exp_takanen2013fig8bookcochleadata.mat');
        for ind=1:length(data.tests)
            data=safe_load('exp_takanen2013fig8bookcochleadata.mat');
            tit=data.tests(ind).case;
            insig=data.tests(ind).cochlear;
            clear data
            % compute the binaural activity map with the model
            output = takanen2013(insig,fs,compType,printFigs,printMap);
            nXBins= length(output.levels)*(size(output.colorMtrx,1)-1);
            dim = size(output.activityMap);
            output.colorGains(output.colorGains>1) =1;
            outputMtrx = single(zeros(dim(1),nXBins,3));
            for colorInd=1:size(output.colorMtrx,1)
                temp = find((output.activityMap==(colorInd-1))==1);
                outputMtrx(temp) = output.colorGains(temp)*output.colorMtrx(colorInd,1);
                outputMtrx(temp+dim(1)*nXBins) = output.colorGains(temp)*output.colorMtrx(colorInd,2);
                outputMtrx(temp+2*dim(1)*nXBins) = output.colorGains(temp)*output.colorMtrx(colorInd,3);
            end
            g(ind)= subplot(3,2,ind);imagesc(output.levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
            clear output outputMtrx
            title(tit);
            set(gca,'YTick',.0:.5:2.5);
            set(gca,'YTickLabel',2.5:-0.5:0);
            set(gca,'Xtick',-1:0.4:1);
            ylabel('Time [s]');
            xlabel('Activation location');
        end
    end
end
%% Figure 9 from the book chapter
if flags.do_fig9
    probDist = zeros(6,19);
    % if the user wishes to compute the cochlear model outputs, binaural
    % input signals are used
    if flags.do_binsig
        data=safe_load('exp_takanen2013fig9bookbinsignals.mat');
        for ind=1:length(data.tests)
            data=safe_load('exp_takanen2013fig9bookbinsignals.mat');
            insig=data.tests(ind).insig;
            tit=data.tests(ind).case;
            clear data
            % compute the binaural activity map with the model
            output = takanen2013(insig,fs,compType,printFigs,printMap);
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
            clear output outputMtrx
            title(tit);
            set(gca,'YTick',1:6);
            set(gca,'YTickLabel',6:-1:1);
            set(gca,'Xtick',-1:0.4:1);
            ylabel('Frequency area');
            xlabel('Distribution of activation');
        end
    end
    %otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea
        data=safe_load('exp_takanen2013fig9bookcochleadata.mat');
        for ind=1:length(data.tests)
            data=safe_load('exp_takanen2013fig9bookcochleadata.mat');
            insig=data.tests(ind).cochlear;
            tit=data.tests(ind).case;
            clear data
            % compute the binaural activity map with the model
            output = takanen2013(insig,fs,compType,printFigs,printMap);
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
            clear output outputMtrx
            title(tit);
            set(gca,'YTick',1:6);
            set(gca,'YTickLabel',6:-1:1);
            set(gca,'Xtick',-1:0.4:1);
            ylabel('Frequency area');
            xlabel('Distribution of activation');
        end
    end
end
%% Figure 7 from takanen2014
if flags.do_fig7_takanen2014
    % compute the cochlear model outputs, load the binaural input signals
    if flags.do_binsig, s='exp_takanen2013fig6artbinsignals.mat'; end
    % otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea, s='exp_takanen2013fig6artcochleadata.mat'; end

    data=safe_load(s);
    data_tests=length(data.tests);
    siglen=zeros(length(data.tests),1);
    data_tests_Data=zeros(length(data_tests),1);
    for ind=1:data_tests
      if flags.do_cochlea
        data_tests_Data(ind)=length(data.tests(ind).cochlearData);
        for caseInd=1:data_tests_Data(ind)
          siglen(ind)=siglen(ind)+length(data.tests(ind).cochlearData(caseInd).cochlear.velocityLeft);
        end
      end
      if flags.do_binsig
        data_tests_Data(ind)=length(data.tests(ind).binSignals);
        for caseInd=1:data_tests_Data(ind)
          siglen(ind)=siglen(ind)+length(data.tests(ind).binSignals(caseInd).insig);
        end
      end          
    end
    clear data % release unused memory
    for ind=1:data_tests
        activityMap=zeros(siglen(ind),114);
        gains=zeros(siglen(ind),114);
        idx=1;
        %some scenarios consist of multiple test cases that are
        %processed separately
        for caseInd=1:data_tests_Data(ind)
            data=safe_load(s);
            if flags.do_cochlea
              insig=data.tests(ind).cochlearData(caseInd).cochlear;
              len=size(insig.velocityLeft,1);
            end
            if flags.do_binsig
              insig=data.tests(ind).binSignals(caseInd).insig;
              len=size(insig,1);
            end            
            ylab=data.tests(ind).ylab;
            scenario=data.tests(ind).scenario;
            ytickPos=data.tests(ind).ytickPos;
            ytickLab=data.tests(ind).ytickLab(end:-1:1);
            clear data % release unused memory
            % compute the binaural activity map with the model
            output = takanen2013(insig,fs,compType,printFigs,printMap);
            %concatenate the separate activity maps into one map
            activityMap(idx:idx+len-1,:)=output.activityMap;
            gains(idx:idx+len-1,:)=output.colorGains;
            idx=idx+len;
            colorMtrx=output.colorMtrx;
            levels=output.levels;
            clear output % release unused memory
        end
        %the anti-phasic sweep contains also frequencies below the
        %frequency range of the model. Hence, the first 0.5 s of the
        %activity map are removed
        if(strcmp('Anti-phasic sinusoidal sweep',scenario)==1)
            activityMap = activityMap(0.5*fs+1:end,:);
            gains = gains(0.5*fs+1:end,:);
        end
        nXBins= length(levels)*(size(colorMtrx,1)-1);
        dim = size(activityMap);
        gains(gains>1) =1;
        outputMtrx = zeros(dim(1),nXBins,3);
        for colorInd=1:size(colorMtrx,1)
            temp = find((activityMap==(colorInd-1))==1);
            outputMtrx(temp) = gains(temp)*colorMtrx(colorInd,1);
            outputMtrx(temp+dim(1)*nXBins) = gains(temp)*colorMtrx(colorInd,2);
            outputMtrx(temp+2*dim(1)*nXBins) = gains(temp)*colorMtrx(colorInd,3);
        end
        g(ind)= subplot(2,2,ind);imagesc(levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
        title(scenario);
        set(gca,'YTick',ytickPos);
        set(gca,'YTickLabel',ytickLab);
        set(gca,'Xtick',-1:0.4:1);
        ylabel(ylab);
        xlabel('Activation location');
    end
end
%% Figure 8 from takanen2014
if flags.do_fig8_takanen2014
    % compute the cochlear model outputs, load the binaural input signals
    if flags.do_binsig, s='exp_takanen2013fig7artbinsignals.mat'; end
    % otherwise pre-computed cochlea model outputs are used
    if flags.do_cochlea, s='exp_takanen2013fig7artcochleadata.mat'; end

    data=safe_load(s);
    data_tests=length(data.tests);
    siglen=zeros(length(data.tests),1);
    data_tests_Data=zeros(length(data_tests),1);
    for ind=1:data_tests
      if flags.do_cochlea
        data_tests_Data(ind)=length(data.tests(ind).cochlearData);
        for caseInd=1:data_tests_Data(ind)
          siglen(ind)=siglen(ind)+length(data.tests(ind).cochlearData(caseInd).cochlear.velocityLeft);
        end
      end
      if flags.do_binsig
        data_tests_Data(ind)=length(data.tests(ind).binSignals);
        for caseInd=1:data_tests_Data(ind)
          siglen(ind)=siglen(ind)+length(data.tests(ind).binSignals(caseInd).insig);
        end
      end          
    end
    clear data % release unused memory
    for ind=1:data_tests
        activityMap=zeros(siglen(ind),114);
        gains=zeros(siglen(ind),114);
        idx=1;
        %some scenarios consist of multiple test cases that are
        %processed separately
        for caseInd=1:data_tests_Data(ind)
            data=safe_load(s);
            if flags.do_cochlea
              insig=data.tests(ind).cochlearData(caseInd).cochlear;
              len=size(insig.velocityLeft,1);
            end
            if flags.do_binsig
              insig=data.tests(ind).binSignals(caseInd).insig;
              len=size(insig,1);
            end            
            ylab=data.tests(ind).ylab;
            scenario=data.tests(ind).scenario;
            ytickPos=data.tests(ind).ytickPos;
            ytickLab=data.tests(ind).ytickLab(end:-1:1);
            clear data % release unused memory
            % compute the binaural activity map with the model
            output = takanen2013(insig,fs,compType,printFigs,printMap);
            %concatenate the separate activity maps into one map
            activityMap(idx:idx+len-1,:)=output.activityMap;
            gains(idx:idx+len-1,:)=output.colorGains;
            idx=idx+len;
            colorMtrx=output.colorMtrx;
            levels=output.levels;
            clear output % release unused memory
        end
        %in order to better visualize the clicks in the precedence
        %effect scenario, most of the silent parts of the signal
        %are removed
        if(strcmp('Precedence effect',scenario)==1)
            activityMap = activityMap([1500:3700 4500:6700 7500:9700 10500:12700 13500:15700 16500:18700 20200:22400],:);
            gains = gains([1500:3700 4500:6700 7500:9700 10500:12700 13500:15700 16500:18700 20200:22400],:);
            gains = 2*gains;
        end
        nXBins= length(levels)*(size(colorMtrx,1)-1);
        dim = size(activityMap);
        gains(gains>1) =1;
        outputMtrx = zeros(dim(1),nXBins,3);
        for colorInd=1:size(colorMtrx,1)
            temp = find((activityMap==(colorInd-1))==1);
            outputMtrx(temp) = gains(temp)*colorMtrx(colorInd,1);
            outputMtrx(temp+dim(1)*nXBins) = gains(temp)*colorMtrx(colorInd,2);
            outputMtrx(temp+2*dim(1)*nXBins) = gains(temp)*colorMtrx(colorInd,3);
        end
        g(ind)= subplot(2,2,ind);imagesc(levels./90,((dim(1)-1):-20:0)/fs,outputMtrx(1:20:end,:,:));
        title(scenario);
        set(gca,'YTick',ytickPos);
        set(gca,'YTickLabel',ytickLab);
        set(gca,'Xtick',-1:0.4:1);
        ylabel(ylab);
        xlabel('Activation location');
    end
end
output = g;

function data=safe_load(filename)
  try
      data=load([amtbasepath,'signals',filesep,filename]);
  catch exception
    disp('=============================================================');
    disp('Please load the necessary mat-files from the companying page:');
    disp('   www.acoustics.hut.fi/publications/papers/AMTool2013-bam/  ');
    disp('and place them in the "signals" directory                ');
    disp('=============================================================');
            
      error('Error: mat-file %s not found',filename);
  end