function output=exp_breebaart2001(varargin)
%EXP_BREEBAART2001   Figures from Breebaart et al.(2001a and 2001b)
%   Usage: output = exp_breebaart2001(flag)
%
%   `exp_breebaart2001(flags,... )` reproduces experiments from the 
%         Breebaart et al (2001a and 2001b) papers.
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'afig2'    Reproduce Fig. 2 from Breebaart et al. (2001a)
%
%
%   See also: data_breebaart2001 and data_vandepar1999
%
%   Examples:
%   ---------
%
%   To display Figure 2 of the 2001a paper use :::
%
%     exp_breebaart2001('afig2');
%

%
%   References:
  
%  AUTHOR: Martina Kreuzbichler


warning off;
%% ------ Check input options --------------------------------------------

  definput.import={'amtcache'};
  definput.flags.type = {'missingflag','afig2','afig6','bfig3','bfig6',...
      'fig1_N0S0_vandepar1999'};
  definput.flags.plot = {'plot','noplot'};

  % Parse input options
  [flags,~]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% Breebaart et al. (2001a) fig. 2
if flags.do_afig2
    
    testsigadlo = amtcache('get','afig2',flags.cachemode);
    
    if isempty(testsigadlo)
        definput.import = {'auditoryfilterbank','ihcenvelope','adaptloop','eicell'};
        definput.importdefaults={'fhigh',8000,'ihc_breebaart','adt_breebaart'};

        [flags,keyvals,~,~,~]  = ltfatarghelper({'flow', 'fhigh', ...
                            'basef'},definput,varargin);


        testsig(:,1) = sin(2*pi*(0:(0.1*48000)-1)'*500/48000);
        testsig(:,2) = sin(2*pi*(0:(0.1*48000)-1)'*4000/48000);

        testsig = setdbspl(testsig,70);
        testsig = [testsig; zeros(0.1*48000,2)];

        for counter = 1:2

            testsigoutmiddle(:,counter) = outmiddleartransfunct(testsig(:,counter),48000);

            [testsigaudfilt(:,:,counter), ~] = auditoryfilterbank(testsigoutmiddle(:,counter),48000,'argimport',flags,keyvals);

            testsigihc(:,:,counter) = ihcenvelope(testsigaudfilt(:,:,counter),48000,'argimport',flags,keyvals);

            testsigadlo(:,:,counter) = adaptloop(testsigihc(:,:,counter),48000,'argimport',flags,keyvals);
        end
        
        amtcache('set','afig2',testsigadlo)
        
    end
    
    output = testsigadlo;
    
% elseif flags.do_afig6
%     definput.import = {'auditoryfilterbank','ihcenvelope','adaptloop','eicell'};
%     definput.importdefaults={'fhigh',8000,'ihc_breebaart','adt_breebaart'};
% 
%     [flags,keyvals,flow,fhigh,basef]  = ltfatarghelper({'flow', 'fhigh', ...
%                         'basef'},definput,varargin);
%                     
%     testsig(:,1) = bandpassnoisefreq(2000,4800, par.ndur,par.nl,par.nbw);
%     
%     testsigoutmiddle(:,counter) = outmiddleartransfunct(testsig(:,counter),48000);
%     
%     [testsigaudfilt(:,:,counter), ~] = auditoryfilterbank(testsigoutmiddle(:,counter),48000,'argimport',flags,keyvals);
%     
%     testsigihc(:,:,counter) = ihcenvelope(testsigaudfilt(:,:,counter),48000,'argimport',flags,keyvals);
% 
%     testsigadlo(:,:,counter) = adaptloop(testsigihc(:,:,counter),48000,'argimport',flags,keyvals);
% 
%     %[siglen,nfreqchannels,naudiochannels,nsignals] = size(testsigadlo);
%     %ei_map = zeros(siglen, nfreqchannels, nsignals);
% 
%     taucount = 1;
%     ildcount = 1;
%     for tau = -0.002:0.0001:0.002
%         for ild = -10:0.5:10
%             ei_map_(taucount,ildcount,:) = eicell(squeeze(testsigadlo(:,9,:,1)),fs,tau,ild); % 9 ist Filter bei 500 Hz
%             ildcount = ildcount + 1;
%         end
%         ildcount = 1;
%         taucount = taucount + 1;
%     end

elseif flags.do_bfig3
    
    [N0Spi125,N0Spi250,N0Spi500,N0Spi1000,N0Spi2000,N0Spi4000] = ...
        amtcache('get','bfig3',flags.cachemode);
    
    if isempty(N0Spi125)
    % do computation
        
        parout = [];
        expset = {'intnum',3,'rule',[2 1],'stepin',8,'stepr',[0.5 2],...
            'stepmin',[1 8],'expstart',65};
        parout = amtafcexp('expinit',parout,expset);

        modelset = {'name','breebaart2001preprocbimono','input1','expsignal',...
            'input2',32000,'input3',0,'input4',0,'outputs',[1 2 3]};
        parout = amtafcexp('modelinit',parout,modelset);

        decisionset = {'name','breebaart2001centralproc','input1','modelout',...
            'input2','modelout', 'input3','modelout'};
        parout = amtafcexp('decisioninit',parout,decisionset);

        centerfreq = [125 250 500 1000 2000 4000];
        bw = [5 10 25 50 100 250];
        bwadd = [500 1000 2000 4000];
        nl = 65;

        % loop for all center freq.
        for freqcount = 1:length(centerfreq)

            %loop for all bandwidths
            for bwcount = 1:length(bw)


                signalset = {'name','breebaart2001siggen','input1','inttyp',...
                    'input2',centerfreq(freqcount),'input3','expvar',...
                    'input4',0.3,'input5',pi,'input6',...
                    bw(bwcount),'input7',nl,'input8',0.4,'input9',0,...
                    'input10',0.05,'input11', 32000};
                parout = amtafcexp('signalinit',parout,signalset);

                % loop for experimental runs
                for runcounter = 1:6
                    result = amtafcexp('run',parout);
                    resultbwvec(runcounter) = result(1)-nl;
                end

                resultvec(bwcount) = mean(resultbwvec);
                resultvecstd(bwcount) = std(resultbwvec,1);
                resultbwvec = zeros(6,1);
                waitbar(bwcount/length(bw)) 

            end

            N0Spi_temp =sprintf('N0Spi%i',centerfreq(freqcount));
            output.(N0Spi_temp) = [resultvec, resultvecstd];

            if freqcount <= length(bwadd)
                bw(end+1) = bwadd(bwcount);
            end

        end
    end
    
    
    


elseif flags.do_fig1_N0S0_vandepar1999
    
    [N0S0125,N0S0250,N0S0500,N0S01000,N0S02000,N0S04000] = ...
        amtcache('get','fig1_N0S0_vandepar1999',flags.cachemode);
    
    if isempty(N0S0125)
%         RUN EXPERIMENT
    end

    output = struct('N0S0125',N0S0125,'N0S0250',N0S0250,'N0S0500',N0S0500,...
        'N0S01000',N0S01000,'N0S02000',N0S02000,'N0S040000',N0S04000);  
end


if flags.do_plot
    if flags.do_afig2
        fax = 0:0.2/9599:0.2;
        hafig2 = figure;
        set(hafig2, 'Position', [400 400 1200 450])
        subplot(1,2,1)
        plot(fax,output(:,9,1))
        axis([0 0.2 -1000 8500])
        xlabel('time [s]')
        ylabel('Output [MU')
        text('Position',[0.1 7000],'string','500 Hz tone')
        subplot(1,2,2)
        plot(fax,output(:,25,2))
        axis([0 0.2 -1000 8500])
        xlabel('time [s]')
        ylabel('Output [MU]')
        text('Position',[0.1 7000],'string','4000 Hz tone')
        suptitle('Fig. 2 Output of the peripheral preprocessor')
        
    elseif flags.do_fig1_N0S0_vandepar1999
        
        %get experimental data
        N0S0expdata125 = data_vandepar1999('fig1_N0S0','nfc125');
        N0S0expdata250 = data_vandepar1999('fig1_N0S0','nfc250');
        N0S0expdata500 = data_vandepar1999('fig1_N0S0','nfc500');
        N0S0expdata1000 = data_vandepar1999('fig1_N0S0','nfc1000');
        N0S0expdata2000 = data_vandepar1999('fig1_N0S0','nfc2000');
        N0S0expdata4000 = data_vandepar1999('fig1_N0S0','nfc4000');
        
        hfig1=figure;
        set(hfig1, 'Position', [100 100 1100 700])
        
        bw = [5 10 25 50 100 250];
        s1 = subplot(3,2,1);
        p1 = get(s1,'pos');
        errorbar(bw,N0S0125(:,1),N0S0125(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([4 5000 -30 10])
        hold on;
        plot(bw,N0S0expdata125,'-og')
        text('Position',[1000 5],'string','125 Hz')
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500];
        s2 = subplot(3,2,2);
        p2 = get(s2,'pos');
        p2(1) = p1(1) + p1(3);
        set(s2, 'pos', p2);
        errorbar(bw,N0S0250(:,1),N0S0250(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([4 5000 -30 10])
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        hold on;
        plot(bw,N0S0expdata250,'-og')
        text('Position',[1000 5],'string','250 Hz')

        bw = [5 10 25 50 100 250 500 1000];
        s3 = subplot(3,2,3);
        p3 = get(s3,'pos');
        p3(2) = p1(2) - p1(4) + 0.005;
        set(s3, 'pos', p3);
        errorbar(bw,N0S0500(:,1),N0S0500(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([1 10000 -30 10])
        set(gca,'XTick',[]);
        set(gca,'YTick',[-30 -20 -10 0]);
        hold on;
        plot(bw,N0S0expdata500,'-og')
        text('Position',[1000 5],'string','500 Hz')
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500 1000 2000];
        s4 = subplot(3,2,4);
        p4 = get(s4,'pos');
        p4(1) = p3(1) + p3(3);
        p4(2) = p3(2);
        set(s4, 'pos', p4);
        errorbar(bw,N0S01000(:,1),N0S01000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        axis([1 10000 -30 10])
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        hold on;
        plot(bw,N0S0expdata1000,'-og')
        text('Position',[1000 5],'string','1000 Hz')
        xlabel('Bandwidth [Hz]');

        bw = [5 10 25 50 100 250 500 1000 2000 4000];
        s5 = subplot(3,2,5);
        p5 = get(s5,'pos');
        p5(2) = p3(2) - p3(4) + 0.005;
        set(s5, 'pos', p5);
        errorbar(bw,N0S02000(:,1),N0S02000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[1 10 100 1000 10000]);
        axis([1 10000 -30 10])
        set(gca,'YTick',[-30 -20 -10 0]);
        hold on;
        plot(bw,N0S0expdata2000,'-og')
        text('Position',[1000 5],'string','2000 Hz')
        xlabel('Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        s6 = subplot(3,2,6);
        p6 = get(s6,'pos');
        p6(1) = p5(1) + p5(3);
        p6(2) = p5(2);
        set(s6, 'pos', p6);
        errorbar(bw,N0S04000(:,1),N0S04000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[ 10 100 1000 10000]);
        axis([1 10000 -30 10])
        set(gca,'YTick',[]);
        hold on;
        plot(bw,N0S0expdata4000,'-og')
        text('Position',[1000 5],'string','4000 Hz')
        xlabel('Bandwidth [Hz]');
        
        legend({'modeled data SSL = 90 dB, NL 70 dB  - monaural factor = 0.0003',...
            'experimental data from van de Paar and Kohlrausch (1999)'},...
            'Position' ,[p2(1)-p2(3)*0.2,0.15,0.1,0.05]);
        
        suptitle('Fig. 1 N_0S_0 Thresholds')
        
    elseif flags.do_bfig3
        
        % get experimental data
        N0Spiexpdata125 = data_vandepar1999('fig1_N0Spi','nfc125');
        N0Spiexpdata250 = data_vandepar1999('fig1_N0Spi','nfc250');
        N0Spiexpdata500 = data_vandepar1999('fig1_N0Spi','nfc500');
        N0Spiexpdata1000 = data_vandepar1999('fig1_N0Spi','nfc1000');
        N0Spiexpdata2000 = data_vandepar1999('fig1_N0Spi','nfc2000');
        N0Spiexpdata4000 = data_vandepar1999('fig1_N0Spi','nfc4000');
        
        % get model data
        N0Spimodeldata125 = data_breebaart2001('fig3','nfc125');
        N0Spimodeldata250 = data_breebaart2001('fig3','nfc250');
        N0Spimodeldata500 = data_breebaart2001('fig3','nfc500');
        N0Spimodeldata1000 = data_breebaart2001('fig3','nfc1000');
        N0Spimodeldata2000 = data_breebaart2001('fig3','nfc2000');
        N0Spimodeldata4000 = data_breebaart2001('fig3','nfc4000');

        hfig2=figure;
        set(hfig2, 'Position', [100 100 1100 700])
        
        bw = [5 10 25 50 100 250];
        subplot(3,2,1);
        errorbar(bw,N0Spi125(:,1),N0Spi125(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata125,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        stem(bw,N0Spiexpdata125,'sg','LineStyle','none','MarkerSize',10)
        text('Position',[1000 -10],'string','125 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500];
        subplot(3,2,2);
        errorbar(bw,N0Spi250(:,1),N0Spi250(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata250,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        stem(bw,N0Spiexpdata250,'sg','LineStyle','none')
        text('Position',[1000 -10],'string','250 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500 1000];
        subplot(3,2,3);
        errorbar(bw,N0Spi500(:,1),N0Spi500(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata500,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        stem(bw,N0Spiexpdata500,'sg','LineStyle','none','MarkerSize',10)
        text('Position',[1000 -10],'string','500 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        bw = [5 10 25 50 100 250 500 1000 2000];
        subplot(3,2,4);
        errorbar(bw,N0Spi1000(:,1),N0Spi1000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata1000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        stem(bw,N0Spiexpdata1000,'sg','LineStyle','none','MarkerSize',10)
        text('Position',[1000 -10],'string','1000 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        
        bw = [5 10 25 50 100 250 500 1000 2000 4000];
        subplot(3,2,5);
        errorbar(bw,N0Spi2000(:,1),N0Spi2000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 -5])
        hold on;
        plot(bw,N0Spimodeldata2000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        stem(bw,N0Spiexpdata2000,'sg','LineStyle','none','MarkerSize',10)
        text('Position',[1000 -10],'string','2000 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');
        
        subplot(3,2,6);
        errorbar(bw,N0Spi4000(:,1),N0Spi4000(:,2),'-x','MarkerSize',10);
        set(gca,'xscal','log')
        set(gca,'XTick',[5 10 25 50 100 250 1000 4000]);
        axis([4 5000 -35 5])
        hold on;
        plot(bw,N0Spimodeldata4000,'-sr','MarkerSize',10,'MarkerFaceColor','r');
        plot(bw,N0Spiexpdata4000,'sg','LineStyle','none','MarkerSize',10)
        text('Position',[1000 0],'string','2000 Hz')
        xlabel('Masker Bandwidth [Hz]');
        ylabel('Threshold S/N [dB]');

        suptitle('Fig. 3 N0Spi Thresholds')
        legend({'modeled data SSL = 65 dB, NL 65 dB  - monaural factor = 0.0003',...
            'model data from Breebaart (2001)',...
            'experimental data from van de Paar and Kohlrausch (1999)'},...
            'Position' ,[0.75,0.9,0.1,0.05]);
            end
end
    
    