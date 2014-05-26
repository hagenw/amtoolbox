% DEMO_GAMMATONE Demo for gammatone.m 
%
%   `demo_gammatone(flag)` demonstrates arrays of impulse responses for
%   implementations of gammatone.m. 
%

%% Global parameters 
    
    flow = 100;
    fhigh = 4000;
    fs = 48000;
    fc = erbspacebw(flow,fhigh);
    nchannels = length(fc); 
    betamul=1;
    n=4;
    insig = [1, zeros(1,8191)];
    treal = (1:length(insig))/fs*1000;
    
%     % ERB spaced channels and corresponding center frequencies
%     for ii = 1:length(fc)
%         fprintf(1, 'Channel: %2d   Center frequency: %f  \n', ii, fc(ii))  
%     end;
    
%% All possible arrays of gammatone impulse responses
 
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
    
% %---- allpole real ----  !!! Not scaled correctly !!!  
%     [b,a] = gammatone(fc,fs,'allpole','causalphase','real');
%     outsig_ap_cph_r = 2*real(ufilterbankz(b,a,insig));
%     outsig_ap_cph_r = permute(outsig_ap_cph_r,[3 2 1]);
%  
%     [b,a] = gammatone(fc,fs,'allpole','peakphase','real');
%     outsig_ap_pph_r = 2*real(ufilterbankz(b,a,insig));
%     outsig_ap_pph_r = permute(outsig_ap_pph_r,[3 2 1]);

%---- allpole complex ----
    [b,a] = gammatone(fc,fs,'allpole','causalphase','complex');
    outsig_ap_cph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_cph_c = permute(outsig_ap_cph_c,[3 2 1]);

    [b,a] = gammatone(fc,fs,'allpole','peakphase','complex');
    outsig_ap_pph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_pph_c = permute(outsig_ap_pph_c,[3 2 1]);
    

%% For plot arrangement

    outsig1 = outsig_cl_cph_r;
    type1   = 'classic / causalphased / real';
    outsig2 = outsig_cl_pph_r;
    type2   = 'classic / peakphased / real';
    outsig3 = outsig_cl_cph_c;
    type3   = 'classic / causalphased / complex';
    outsig4 = outsig_cl_pph_c;
    type4   = 'classic / peakphased / complex';
    outsig5 = outsig_ap_cph_c;
    type5   = 'allpole / causalphased / complex';
    outsig6 = outsig_ap_pph_c;
    type6   = 'allpole / peakphased / complex';
       
        
%% Plot;
    
    % figure 1 - classic real;
    figure('units','normalized','outerposition',[0 0.5 1 0.5])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*outsig1(ii,:) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 25]);
    set(gca,'XTick',xt, 'YTick', yt)
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*outsig2(ii,:) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type2)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 25]);
    set(gca,'XTick',xt, 'YTick', yt)
    hold off
    
    % figure 2 - classic complex;
    figure ('units','normalized','outerposition',[0 0.05 1 0.5])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*outsig3(ii,:) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type3)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 25]);
    set(gca,'XTick',xt, 'YTick', yt)
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*outsig4(ii,:) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type4)])        
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 25]);
    set(gca,'XTick',xt, 'YTick', yt)
    hold off
    
    % figure 3 - allpass complex;
    figure ('units','normalized','outerposition',[0 0 1 0.5])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*outsig5(ii,:) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type5)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 25]);
    set(gca,'XTick',xt, 'YTick', yt)
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*outsig6(ii,:) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gamma impulse responses of type: ', num2str(type6)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 25]);
    set(gca,'XTick',xt, 'YTick', yt)
    hold off
    
    
%% Plot backup; 
% 
%     figure('units','normalized','outerposition',[0 0 1 1])
%     subplot(3,2,1)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,25*outsig1(ii,:) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type1)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 25]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     hold off
%         
%     subplot(3,2,2)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,25*outsig2(ii,:) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type2)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 25]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     hold off
%         
%     subplot(3,2,3)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,25*outsig3(ii,:) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type3)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 25]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     hold off
%         
%     subplot(3,2,4)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,25*outsig4(ii,:) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type4)])        
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 25]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     hold off
%         
%     subplot(3,2,5)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,25*outsig5(ii,:) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type5)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])    xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 25]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     hold off
%         
%     subplot(3,2,6)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,25*outsig6(ii,:) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gamma impulse responses of type: ', num2str(type6)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])    xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 25]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     hold off
       