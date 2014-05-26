% DEMO_GAMMATONE Demo for gammatone.m 
%
%   `demo_gammatone(flag)` demonstrates arrays of impulse responses for the
%   gammatone filterbank implementations of gammatone.m. and Hohmanns
%   allpass implementation marked by gfb_.
%
%   Figure 1 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) derived from 
%   gammatone.m with real-valued filter coefficients and in the second plot
%   this implementation with option peakphased.
%
%   Figure 2 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) derived from
%   gammatone.m with complex-valued filter coefficients and in the second
%   plot this implementation with option peakphased.
%
%   Figure 3 shows in the first plot an array of 24 erb-spaced channels of
%   allpole gammatone implementation (Hohmann, 2002) derived from gammatone.m
%   with real-valued filter coefficients and in the second plot this
%   implementation with option peakphased.
%
%   Figure 4 shows in the second plot an array of 24 erb-spaced channels of
%   allpass gammatone implementation (Hohmann, 2002) with complex-valued 
%   filter coefficients.
%
%   See also: exp_gammatone gammatone exp_hohmann2002 demo_hohmann2002
%   gfb_analyzer gfb_analyzer_process

% AUTHOR: CK, 2014

%% Gammtone implementation;

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 48000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    betamul=1;                      % ERB density multiplicator;
    n=4;                            % Filter order;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    
%% All possible arrays of gammatone impulse responses;
 
%---- classic real ----
    % Derive filter coefficients and filter impulse input signal; 
    [b,a] = gammatone(fc,fs,'classic','causalphase','real');
    outsig_cl_cph_r = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_cph_r = permute(outsig_cl_cph_r,[3 2 1]);
    
    % Derive filter coefficients and filter impulse input signal with option peakphase;
    [b,a] = gammatone(fc,fs,'classic','peakphase','real');
    outsig_cl_pph_r = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_pph_r = permute(outsig_cl_pph_r,[3 2 1]);

%---- classic complex ----
    disp('Classic, complex-valued, causal-phased gammtone implementation:')
    % Derive filter coefficients and filter impulse input signal;
    [b,a] = gammatone(fc,fs,'classic','causalphase','complex');
    outsig_cl_cph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_cph_c = permute(outsig_cl_cph_c,[3 2 1]);

    disp('Classic, complex-valued, peak-phased gammtone implementation:')
    % Derive filter coefficients and filter impulse input signal with option peakphase;
    [b,a] = gammatone(fc,fs,'classic','peakphase','complex');
    outsig_cl_pph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_cl_pph_c = permute(outsig_cl_pph_c,[3 2 1]);
    
% %---- allpole real ----  !!! Not scaled correctly !!!
%       % Derive filter coefficients and filter impulse input signal;
%     [b,a] = gammatone(fc,fs,'allpole','causalphase','real');
%     outsig_ap_cph_r = 2*real(ufilterbankz(b,a,insig));
%     outsig_ap_cph_r = permute(outsig_ap_cph_r,[3 2 1]);
%  
%       % Derive filter coefficients and filter impulse input signal with option peakphase;
%     [b,a] = gammatone(fc,fs,'allpole','peakphase','real');
%     outsig_ap_pph_r = 2*real(ufilterbankz(b,a,insig));
%     outsig_ap_pph_r = permute(outsig_ap_pph_r,[3 2 1]);

%---- allpole complex ----
    % Derive filter coefficients and filter impulse input signal;
    [b,a] = gammatone(fc,fs,'allpole','causalphase','complex');
    outsig_ap_cph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_cph_c = permute(outsig_ap_cph_c,[3 2 1]);

    % Derive filter coefficients and filter impulse input signal with option peakphase;
    [b,a] = gammatone(fc,fs,'allpole','peakphase','complex');
    outsig_ap_pph_c = 2*real(ufilterbankz(b,a,insig));
    outsig_ap_pph_c = permute(outsig_ap_pph_c,[3 2 1]);

    
%% Hohmann implementation;

    % Parameters;
    fs = 48000;                 % Sampling rate in Hz;
    flow = 100;                 % Lowest center frequency in Hz;
    basef = 888.44;             % Base center frequency in Hz;
    fhigh  = 4000;              % Highest center frequency in Hz;
    filters_per_ERBaud = 1;     % Filterband density on ERB scale;     
    gamma_order= 4;             % Filter order;
    bandwidth_factor = 1;       % Bandwidth factor;
    
    % Construct new analyzer object;
    analyzer = gfb_analyzer_new(fs,flow, basef, fhigh,filters_per_ERBaud,gamma_order,bandwidth_factor);
    % Impulse signal;
    impulse = [1, zeros(1,8191)];
    % Filter signal;
    [impulse_response, analyzer] = gfb_analyzer_process(analyzer, impulse);

    
%% For plot arrangement;
    
    % Figure 1;
    outsig1 = outsig_cl_cph_r;                    % Outputsinal 1;
    type1   = 'classic / causalphased / real';    % Type of implemantion for headline;
    outsig2 = outsig_cl_pph_r;                    % Outputsinal 2;
    type2   = 'classic / peakphased / real';      % Type of implemantion for headline;
    
    % Figure 2;
    outsig3 = outsig_cl_cph_c;                    % Outputsinal 1;
    type3   = 'classic / causalphased / complex'; % Type of implemantion for headline;
    outsig4 = outsig_cl_pph_c;                    % Outputsinal 2;
    type4   = 'classic / peakphased / complex';   % Type of implemantion for headline;

    % Figure 3;
    outsig5 = outsig_ap_cph_c;                    % Outputsinal 1;
    type5   = 'allpole / causalphased / complex'; % Type of implemantion for headline;
    outsig6 = outsig_ap_pph_c;                    % Outputsinal 2;
    type6   = 'allpole / peakphased / complex';   % Type of implemantion for headline;
    
    % Figure 4;
    outsig7 = impulse_response;                             % Outputsinal 1;
    type7   = 'hohmann / allpole / peakphased / complex';   % Type of implemantion for headline;
       
    
%% Plot;
    
    % Figure 1 - classic real;
    figure('units','normalized','outerposition',[0 0.525 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*real(outsig1(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(fc(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(fc(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(fc(8),'%6.2f') ' Hz'], ['#12:  ' num2str(fc(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(fc(16),'%6.2f') ' Hz'], ['#20: ' num2str(fc(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(fc(24),'%6.2f') ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*real(outsig2(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type2)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(fc(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(fc(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(fc(8),'%6.2f') ' Hz'], ['#12:  ' num2str(fc(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(fc(16),'%6.2f') ' Hz'], ['#20: ' num2str(fc(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(fc(24),'%6.2f') ' Hz']})
    box on
    hold off
    
    
    % Figure 2 - classic complex;
    figure ('units','normalized','outerposition',[0 0.05 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*real(outsig3(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type3)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(fc(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(fc(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(fc(8),'%6.2f') ' Hz'], ['#12:  ' num2str(fc(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(fc(16),'%6.2f') ' Hz'], ['#20: ' num2str(fc(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(fc(24),'%6.2f') ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*real(outsig4(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type4)])        
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(fc(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(fc(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(fc(8),'%6.2f') ' Hz'], ['#12:  ' num2str(fc(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(fc(16),'%6.2f') ' Hz'], ['#20: ' num2str(fc(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(fc(24),'%6.2f') ' Hz']})
    box on
    hold off
    
    % Figure 3 - allpass complex;
    figure ('units','normalized','outerposition',[0 0.525 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*real(outsig5(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type5)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(fc(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(fc(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(fc(8),'%6.2f') ' Hz'], ['#12:  ' num2str(fc(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(fc(16),'%6.2f') ' Hz'], ['#20: ' num2str(fc(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(fc(24),'%6.2f') ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*real(outsig6(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gamma impulse responses of type: ', num2str(type6)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(fc(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(fc(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(fc(8),'%6.2f') ' Hz'], ['#12:  ' num2str(fc(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(fc(16),'%6.2f') ' Hz'], ['#20: ' num2str(fc(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(fc(24),'%6.2f') ' Hz']})
    box on
    hold off
    
    % Figure 4 - hohmann allpass complex;
    figure ('units','normalized','outerposition',[0.475 0.05 0.475 0.475])
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,30*real(outsig7(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type7)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  ' num2str(analyzer.center_frequencies_hz(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(analyzer.center_frequencies_hz(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(analyzer.center_frequencies_hz(8),'%6.2f') ' Hz'], ['#12:  ' num2str(analyzer.center_frequencies_hz(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(analyzer.center_frequencies_hz(16),'%6.2f') ' Hz'], ['#20: ' num2str(analyzer.center_frequencies_hz(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(analyzer.center_frequencies_hz(24),'%6.2f') ' Hz']})
    box on
    hold off
        
    