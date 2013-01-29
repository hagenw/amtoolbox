function exp_dietz2011(varargin)
%EXP_DIETZ2011 Experiments from Dietz et al. 2011
%
%   `exp_dietz2011(fig)` reproduce Fig. no. *fig* from the Dietz et
%   al. 2011 paper.
%
%   **Note**: The input signals used in this routine are not identical to
%   the ones used in the original paper.
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'fig3'    Reproduce Fig. 3 panels a + b.
%
%     'fig4'    Reproduce Fig. 4.
%
%     'fig5'    Reproduce Fig. 5.
%
%     'fig6'    Reproduce Fig. 6.
%
%   Examples:
%   ---------
%
%   To display Fig. 3 use :::
%
%     exp_dietz2011('fig3');
%
%   To display Fig. 4 use :::
%
%     exp_dietz2011('fig4');
%
%   To display Fig. 5 use :::
%
%     exp_dietz2011('fig5');
%
%   To display Fig. 6 use :::
%
%     exp_dietz2011('fig6');
%
%   See also: dietz2011
%
%   References: dietz2011auditory

% AUTHOR: Mathias Dietz

definput.flags.type = {'missingflag','fig3','fig4','fig5','fig6'};
definput.flags.plot = {'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
               sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

if flags.do_fig3

    load azi_lookup_4;
    [signal,fs]=data_dietz2011('five_speakers');
    s_pos = [-80 -30 0 30 80];
    ic_threshold=0.98;
    panellabel = 'ab';

    % run IPD model on signal
    [hairc_fine, hairc_mod, fc, hairc_ild]=dietz2011(signal,fs);

    % convert interaural information into azimuth
    angl=itd2angle(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,lookup,2.5);

    h_ic=zeros(91,12);
    h_all=histc(angl,-90:2:90);
    for n=1:12
        h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-90:2:90);
    end

    if flags.do_plot
        figure;
        fontsize = 14;
        set(gcf,'Position',[100 100 1170 700])
        
        for panel = 1:2
            subplot(1,2,panel)
            switch panel
                
              case 1
                bar(-90:2:90,sum(h_all,2),'r')
                title('Mean histogram of all fine-structure channels','Fontsize',fontsize)
                ymax = max(sum(h_all,2));
              case 2
            bar(-90:2:90,sum(h_ic,2))
            title('Mean histogram with VS filter','Fontsize',fontsize)
            ymax = max(sum(h_ic,2));
            end
            set(gca,'Fontsize',fontsize)
            set(gca,'XTick',s_pos)
            xlim([-93 93])
            ylim([0 ymax*1.1])
            xlabel('Azimuth [deg]','Fontsize',fontsize)
            ylabel('Frequency of occurence','Fontsize',fontsize)
            text (-80,ymax*.95,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
            
        end
    end;
end;


if flags.do_fig4
    
    % This reproduces Figure 4 from Dietz et al Speech Comm. 2011

    load azi_lookup_4;
    signal=data_dietz2011('two_speakers');
    fs = 16000;
    s_pos = [-80 -30 0 30 80];
    ic_threshold=0.98;
    cn = [10 1]; % channel numbers for separate plots (1st entry also for time plot)
    panellabel = 'abc';

    % run IPD model on signal
    disp('mod_center_frequency_hz could have been set to either 135 or 216')
    [hairc_fine, hairc_mod, fc, hairc_ild]=dietz2011(signal,fs);
    % convert interaural information into azimuth
    angl=itd2angle(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,lookup,2.5);
    angl_fmod=hairc_mod.itd_lp(:,13:23)*140000; %linear approximation. paper version is better than this

    h_ic=zeros(61,12);
    h_all=histc(angl,-60:2:60);
    h_fmod=histc(nonzeros(angl_fmod),-60:2:60);
    for n=1:12
        h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-60:2:60);
    end

    if flags.do_plot
        figure;
        fontsize = 14;
        set(gcf,'Position',[100 100 1170 400])
        
        for panel = 1:3
            subplot(1,3,panel)
            switch panel
                
              case 1
                bar(-60:2:60,sum(h_fmod,2))
                title('histogram of mod ITD channels 13-23','Fontsize',fontsize)
                ymax = max(sum(h_fmod,2));
              case 2
                axis;
                text (-60,ymax*.7,'mod center frequency hz has to be','Fontsize',fontsize-3)
            text (-60,ymax*.6,'set manually to either 135 or 216','Fontsize',fontsize-3)
            text (-60,ymax*.5,'result will be displayed in panel a','Fontsize',fontsize-3)
              case 3
                bar(-60:2:60,sum(h_ic,2))
                title('histogram of fine ITD channels 1-12','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
            end
            set(gca,'Fontsize',fontsize)
            set(gca,'XTick',s_pos)
            xlim([-63 63])
            ylim([0 ymax*1.1])
            xlabel('Azimuth [deg]','Fontsize',fontsize)
            ylabel('Frequency of occurence','Fontsize',fontsize)
            text (-50,ymax*.97,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
            
        end
    end;

end;

if flags.do_fig5
    % mix signals
    signal1=data_dietz2011('one_of_three');
    signal2=data_dietz2011('two_of_three');
    signal3=data_dietz2011('three_of_three');
    noise=wavread('bNoise2011');
    noise = noise(1:40000,:);
    fs = 16000;

    % derive histograms
    ic_threshold=0.98;
    load azi_lookup_4;
    k = zeros(14,38);
    for n = 1:7
        switch n
          case 1
            signal = signal1+noise;
            
          case 2
            signal = signal2+noise;
          case 3
            signal = signal3+noise;
          case 4
            signal = signal1+2*noise;
          case 5
            signal = signal2+2*noise;
          case 6
            signal = signal3+2*noise;
          case 7
            signal = noise;
        end
        % run IPD model on signal
        [hairc_fine, hairc_mod, fc, hairc_ild]=dietz2011(signal,fs);
        % convert interaural information into azimuth
        angl=itd2angle(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,lookup,2.5);

        h_ic=zeros(38,12);
        k(2*n,:)=sum(histc(angl,-92.5:5:92.5),2);
        for erb=1:12
            h_ic(:,erb)=histc(angl(hairc_fine.ic(:,erb)>ic_threshold&[diff(hairc_fine.ic(:,erb))>0; 0],erb),-92.5:5:92.5);
        end
        k(2*n-1,:)=sum(h_ic,2);
    end

    if flags.do_plot
        % plot
        figure;
        set(gcf,'Position',[100 100 980 500])
        y=[0.14 0.57];
        x=[0.1 0.22 0.34 0.49 0.61 0.73 0.88];
        cols = 2;
        condi={'1S 0dB','','2S 0dB','',...
               '3S 0dB','','1S -6dB','',...
               '2S -6dB','','3S -6dB','','noise',''};
        
        for n = 1:14
            axes('position',[x(ceil(n/cols)) y(mod(n-1,cols)+1) 0.11 0.42],...
                 'box','on','fontSize',12)
            
            if mod(n,2)==0
                set(gca,'YTick',[10 20 30],'YTickLabel',{'','',''},'ylim',[0 32])
                set(gca,'Visible','on','XTick',[-30 0 30],...
                        'xlim',[-50 50],'XTickLabel',{'','','',''})
            else
                set(gca,'YTick',[1 2 3],'YTickLabel',{'','',''},'ylim',[0 3.33])
                set(gca,'Visible','on','XTick',[-30 0 30],...
                        'xlim',[-50 50],'XTickLabel',{'-30','0','+30'})
                xlabel(condi(n))
            end
            
            if n==2
                set(gca,'YTick',[10 20 30],'YTickLabel',{'10k','20k','30k'})
                ylabel('no IVS mask')
            elseif n==1        
                set(gca,'YTick',[1 2 3],'YTickLabel',{'1k','2k','3k'})
                ylabel('with IVS mask')
            end
            hold on
            bar(-92.5:5:92.5,k(n,:)/1000)
            colormap gray 
            
        end
        hold off
    end;
end;

if flags.do_fig6

    load azi_lookup_4;
    signal=data_dietz2011('one_speaker_reverb');
    fs = 16000;
    s_pos = [-45 0 45];
    ic_threshold=0.98;
    cn = [10 1]; % channel numbers for separate plots (1st entry also for time plot)
    panellabel = 'abcd';

    % run IPD model on signal
    [hairc_fine, hairc_mod, fc, hairc_ild]=dietz2011(signal,fs);
    % convert interaural information into azimuth
    angl=itd2angle(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,lookup,2.5);
    angl_fmod=hairc_mod.itd_lp(:,13:23)*140000; %linear approximation. paper version is better than this

    h_ic=zeros(71,12);
    h_all=histc(angl,-70:2:70);
    h_fmod=histc(nonzeros(angl_fmod),-70:2:70);
    for n=1:12
        h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-70:2:70);
    end

    if flags.do_plot
        figure;
        fontsize = 14;
        set(gcf,'Position',[60 100 1370 350])
        
        for panel = 1:4
            subplot(1,4,panel)
            switch panel
                
              case 1
                bar(-70:2:70,sum(h_ic,2))
                title('fine ITD channels 200-1400 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
              case 2
                bar(-70:2:70,sum(h_ic(:,6:12),2))
                title('fine ITD channels 500-1400 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
              case 3
                bar(-70:2:70,sum(h_ic(:,1:3),2))
                title('fine ITD channels 200-400 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
              case 4
                bar(-70:2:70,sum(h_fmod,2))
                title('mod ITD channels 1400-5000 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
            end
            set(gca,'Fontsize',fontsize)
            set(gca,'XTick',s_pos)
            xlim([-73 73])
            ylim([0 ymax*1.1])
            xlabel('Azimuth [deg]','Fontsize',fontsize)
            ylabel('Frequency of occurence','Fontsize',fontsize)
            text (-50,ymax*.97,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
            
        end        
    end;    
end;