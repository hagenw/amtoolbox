% DEMO_GAMMATONE Demo for gammatone.m and gammatone fiterbank gfb_hohmann
%
%   `demo_gammatone(flag)` demonstrates arrays of impulse responses for all
%   implementations of gammatone.m. plus , compares different implematations for 
%   selected frequency band, and compares complex allpole implementation
%   with Hohmanns gfb_hohmann filterbank.
%
%
%

%% Global parameters
    
    flow = 20;
    fhigh = 8000;
    fs = 24000;
    fc = erbspacebw(flow,fhigh);
    nchannels = length(fc); 
    betamul=1;
    n=4;
    insig = [1, zeros(1,255)];
    
    % ERB spaced channels and corresponding center frequencies
    for ii = 1:length(fc)
        fprintf(1, 'Channel: %2d   Center frequency: %f  \n', ii, fc(ii))  
    end;
    
%% All array figures
 
 %---- classic real ----
    [b,a] = gammatone(fc,fs,'classic','causalphase','real');
    outsig_cl_cph_r = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_cph_r = permute(outsig_cl_cph_r,[3 2 1]);
    
    [b,a] = gammatone(fc,fs,'classic','peakphase','real');
    outsig_cl_pph_r = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_pph_r = permute(outsig_cl_pph_r,[3 2 1]);

%---- classic complex ----    
    [b,a] = gammatone(fc,fs,'classic','causalphase','complex');
    outsig_cl_cph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_cph_c = permute(outsig_cl_cph_c,[3 2 1]);
    
    [b,a] = gammatone(fc,fs,'classic','peakphase','complex');
    outsig_cl_pph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_pph_c = permute(outsig_cl_pph_c,[3 2 1]);
    
%---- allpass real ----    
    [b,a] = gammatone(fc,fs,'allpole','causalphase','real');
    outsig_ap_cph_r = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_cph_r = permute(outsig_ap_cph_r,[3 2 1]);
 
    [b,a] = gammatone(fc,fs,'allpole','peakphase','real');
    outsig_ap_pph_r = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_pph_r = permute(outsig_ap_pph_r,[3 2 1]);

%---- allpass complex ----
    [b,a] = gammatone(fc,fs,'allpole','causalphase','complex');
    outsig_ap_cph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_cph_c = permute(outsig_ap_cph_c,[3 2 1]);
 
    [b,a] = gammatone(fc,fs,'allpole','peakphase','complex');
    outsig_ap_pph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_pph_c = permute(outsig_ap_pph_c,[3 2 1]);
    

%% Classic vs Allpole causalphased


    %if flags.do_classicVsAllpole 
        outsig1 = outsig_cl_cph_r;
        type1   = 'classic / causalphased / real';
        outsig2 = outsig_cl_cph_c;
        type2   = 'classic / causalphased / complex';
        outsig3 = outsig_ap_cph_r;
        type3   = 'allpole / causalphased / real';
        outsig4 = outsig_ap_cph_c;
        type4   = 'allpole / causalphased / complex';
    %end;
    
    figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,2,1)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig1(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type1)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,2)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig2(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type2)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,3)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig3(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type3)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,4)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig4(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type4)])        
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/  ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
%% Classic vs Allpole peakphased

    % if flags.do_classicVsAllpolePeakphased
        outsig1 = outsig_cl_pph_r;
        type1   = 'classic / peakphased / real';
        outsig2 = outsig_cl_pph_c;
        type2   = 'classic / peakphased / complex';
        outsig3 = outsig_ap_pph_r;
        type3   = 'allpole / peakphased / real';
        outsig4 = outsig_ap_pph_c;
        type4   = 'allpole / peakphased / complex';
    % end;
    
    figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,2,1)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig1(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type1)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,2)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig2(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type2)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,3)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig3(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type3)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,4)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig4(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type4)])        
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/  ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off

%% Classic causalphased vs peakphased        
    
    % if flags.do_classicCausalphasedVsPeakphased    
        outsig1 = outsig_cl_cph_r;
        type1   = 'classic / causalphased / real';
        outsig2 = outsig_cl_cph_c;
        type2   = 'classic / causalphased / complex';
        outsig3 = outsig_cl_pph_r;
        type3   = 'classic / peakphased / real';
        outsig4 = outsig_cl_pph_c;
        type4   = 'classic / peakphased / complex';
    % end;
    
    figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,2,1)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig1(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type1)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,2)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig2(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type2)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,3)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig3(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type3)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,4)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig4(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses of type ', num2str(type4)])        
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/  ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off

%% Allpass causalphased vs peakphased

    %if flags.do_allpassCausalphasedVsPeakphased    
        outsig1 = outsig_ap_cph_r;
        type1   = 'allpass / causalphased / real';
        outsig2 = outsig_ap_cph_c;
        type2   = 'allpass / causalphased / complex';
        outsig3 = outsig_ap_pph_r;
        type3   = 'allpass / peakphased / real';
        outsig4 = outsig_ap_pph_c;
        type4   = 'allpass / peakphased / complex';
    %end;
       
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,2,1)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig1(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel gammatone impulse responses of type ', num2str(type1)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,2)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig2(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channels gammatone impulse responses of type ', num2str(type2)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,3)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig3(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channels gammatone impulse responses of type ', num2str(type3)])
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/ ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        subplot(2,2,4)
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(14*outsig4(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channels gammatone impulse responses of type ', num2str(type4)])        
        xlabel 'Time (ms)'
        ylabel ([ num2str(nchannels),' channels/  ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
        %set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        