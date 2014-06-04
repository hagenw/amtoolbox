% DEMO_GAMMATONE Demo for gammatone.m 
%
%   `demo_gammatone(flag)` demonstrates arrays of impulse responses for the
%   gammatone filterbank implementations of gammatone.m (classic and allpass
%   implementations) and Hohmanns allpass implementation marked by gfb_.
%   
%   .. figure::
%
%      Classic gammatone implementation with real-valued filter coefficients derived from gammatone.m.
%
%   Figure 1 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) derived from 
%   gammatone.m with real-valued filter coefficients and in the second plot
%   this implementation with option 'peakphase', which makes the phase of each
%   filter be zero when the envelope of the impulse response of the filter peaks.
%
%   .. figure::
%
%      Classic gammatone implementation with complex-valued filter coefficients derived from gammatone.m.
%
%   Figure 2 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) derived from
%   gammatone.m with complex-valued filter coefficients and in the second
%   plot this implementation with option 'peakphase', which makes the phase
%   of each filter be zero when the envelope of the impulse response of the
%   filter peaks (both do not scale correctly).
%
%   .. figure::
%
%      Allpole gammatone implementation with complex-valued filter coefficients derived from gammatone.m.
%
%   Figure 3 shows in the first plot an array of 24 erb-spaced channels of
%   allpole gammatone implementation (Hohmann, 2002) derived from gammatone.m
%   with complex-valued filter coefficients and in the second plot this
%   implementation with option 'peakphase', which makes the phase of each
%   filter be zero when the envelope of the impulse response of the filter peaks.
%
%   .. figure::
%
%      Allpole gammatone implementation with complex-valued filter coefficients derived from gfb_.
%
%   Figure 4 shows an array of 24 erb-spaced channels of allpass gammatone
%   implementation (Hohmann, 2002) with complex-valued filter coefficients
%   derived from gfb_ and in the second plot this implementation peakphased,
%   which makes the phase of each filter be zero when the envelope of the
%   impulse response of the filter peaks and delayed so the peaks of each
%   filter are arranged above each other.
%
%   .. figure::
%
<<<<<<< HEAD
%      Classic gammatone implementation with real-valued filter coefficients derived from gammatone.m with otpion '6dBperoctave'.
=======
%      Classic gammatone implementation with real-valued filter coefficients
%      derived from gammatone.m with otpion '6dBperoctave'.
>>>>>>> affe3613355ea81f48e6e373b6339197f89b3a80
%
%   Figure 5 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) with real-valued
%   filter coefficients and option '6dBperoctave', which scales the amplitude
%   +/- 6 dB per octave with 0 dB at 4000 Hz and in the second plot this
%   implementation with option 'peakphase_new', which makes the phase of each
%   filter be zero when the envelope of the impulse response of the filter peaks.
%
%   .. figure::
%
<<<<<<< HEAD
%      Allpole gammatone implementation with real-valued filter coefficients derived from gammatone.m with otpion '6dBperoctave'.
=======
%      Allpole gammatone implementation with real-valued filter coefficients
%      derived from gammatone.m with otpion '6dBperoctave'.
>>>>>>> affe3613355ea81f48e6e373b6339197f89b3a80
%
%   Figure 6 shows in the first plot an array of 24 erb-spaced channels of
%   allpole gammatone implementation (Hohmann, 2002) with complex-valued
%   filter coefficients and option '6dBperoctave', which scales the amplitude
%   +/- 6 dB per octave with 0 dB at 4000 Hz and in the second plot this
%   implementation with option 'peakphase_new', which makes the phase of each
%   filter be zero when the envelope of the impulse response of the filter peaks.
%
%   .. figure::
%
<<<<<<< HEAD
%      Classic with real-valued and allpole with complex-valued filter coefficients of gammatone implementation with otpion '6dBperoctave', peakphased and delayed.
=======
%      Classic with real-valued and allpole with complex-valued filter
%      coefficients of gammatone implementation with otpion '6dBperoctave',
%      peakphased and delayed.
>>>>>>> affe3613355ea81f48e6e373b6339197f89b3a80
%
%   Figure 7 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) with real-valued
%   filter coefficients and option '6dBperoctave', which scales the amplitude
%   +/- 6 dB per octave with 0 dB at 4000 Hz and option 'peakphase_new', which
%   makes the phase of each filter be zero when the envelope of the impulse
%   response of the filter peaks. Further the filters are delayed so the peaks
%   of each filter are arranged above each other. The second plot shows the
%   same but with allpole implementation (Hohmann, 2002) with complex-valued
%   filter coefficients derived from gammatone.m.   
%
%   See also: exp_gammatone gammatone exp_hohmann2002 demo_hohmann2002
%   gfb_analyzer_new gfb_analyzer_process

% AUTHOR: Christian Klemenschitz, 2014

%% Figure 1

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    %---- classic real ----
    % Derive filter coefficients and filter impulse responses; 
    [b,a] = gammatone(fc,fs,'classic','causalphase','real');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);
    
    % Derive filter coefficients and filter impulse responses with option 'peakphase';
    [b,a] = gammatone(fc,fs,'classic','peakphase','real');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Figure 1;
    type1   = 'classic / causalphased / real';    % Type of implemantion for headline;
    type2   = 'classic / peakphased / real';      % Type of implemantion for headline;
    
    % Plot
    % Figure 1 - classic real;
    figure('units','normalized','outerposition',[0 0.525 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig1(ii,:)) + dy)
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
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig2(ii,:)) + dy)
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
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
%% Figure 2

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    %---- classic complex ----
    disp('Classic, complex-valued, causal-phased gammtone implementation:')
    % Derive filter coefficients and filter impulse responses;
    [b,a] = gammatone(fc,fs,'classic','causalphase','complex');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);

    disp('Classic, complex-valued, peak-phased gammtone implementation:')
    % Derive filter coefficients and filter impulse responses with option 'peakphase';
    [b,a] = gammatone(fc,fs,'classic','peakphase','complex');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Figure 2;
    type1   = 'classic / causalphased / complex'; % Type of implemantion for headline;
    type2   = 'classic / peakphased / complex';   % Type of implemantion for headline;
    
    % Plot
    % Figure 2 - classic complex;
    figure('units','normalized','outerposition',[0 0.05 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig1(ii,:)) + dy)
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
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig2(ii,:)) + dy)
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
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
%% Figure 3

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    %---- allpole complex ----
    % Derive filter coefficients and filter impulse responses;
    [b,a] = gammatone(fc,fs,'allpole','causalphase','complex');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);

    % Derive filter coefficients and filter impulse responses with option 'peakphase';
    [b,a] = gammatone(fc,fs,'allpole','peakphase','complex');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Figure 3;
    type1   = 'allpole / causalphased / complex'; % Type of implemantion for headline;
    type2   = 'allpole / peakphased / complex';   % Type of implemantion for headline;
    
    % Plot
    % Figure 3 - allpole complex;
    figure('units','normalized','outerposition',[0 0.05 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig1(ii,:)) + dy)
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
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig2(ii,:)) + dy)
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
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
%% Figure 4

    % Hohmann implementation;
    % Parameters;
    fs = 28000;                 % Sampling rate in Hz;
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
       
    % Find peak at envelope maximum for lowest channel and add one sample;
    delay_samples = find(abs(impulse_response(1,:)) == max(abs(impulse_response(1,:)))) + 1;
    
    % 
    delay = gfb_delay_new(analyzer, delay_samples);
     
    [outsig, delay] = gfb_delay_process(delay, impulse_response);
    
    % Figure 4;
    type1   = 'hohmann / allpole / complex';   % Type of implemantion for headline;
    type2   = 'hohmann / allpole / complex / peakphased and delayed';   % Type of implemantion for headline;
    
    % Plot
    % Figure 4 - hohmann allpass complex;
    figure ('units','normalized','outerposition',[0 0.05 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(impulse_response(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig(ii,:)) + dy)
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
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
   
%% Figure 5

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    
    %---- classic real 6dBperoctave ----    
    % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
    [b,a] = gammatone(fc,fs,'classic','real','6dBperoctave');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);
    
    %---- classic real peakphase(new) 6dBperoctave ----
    % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
    [b,a] = gammatone(fc,fs,'classic','real','peakphase_new','6dBperoctave');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Figure 5;
    type1   = 'classic / casualphased / real / +/-6dB';   % Type of implemantion for headline;
    type2   = 'classic / peakphased(new) / real / +/-6dB';  % Type of implemantion for headline;
    
    % Plot
    % Figure 5 - classic / causalphased / real / +/-6dB;
    figure('units','normalized','outerposition',[0 0.525 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,8*real(outsig1(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, -4, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08: ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,8*real(outsig2(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type2)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, -4, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off


%% Figure 6

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
   
    %---- allpole complex 6dBperoctave ----    
    % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
    [b,a] = gammatone(fc,fs,'allpole','complex','6dBperoctave');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);
    
    %---- allpole peakphase(new) complex 6dBperoctave ----
    % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
    [b,a] = gammatone(fc,fs,'allpole','complex','peakphase_new','6dBperoctave');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Figure 6;
    type1   = 'allpole / causalphased  / complex / +/-6dB';   % Type of implemantion for headline;
    type2   = 'allpole / peakphased(new) / complex / +/-6dB';  % Type of implemantion for headline;
    
    % Plot
    % Figure 6 - classic / causalphased / real / +/-6dB;
    figure('units','normalized','outerposition',[0 0.05 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,8*real(outsig1(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, -2, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,8*real(outsig2(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type2)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, -2, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
%% Figure 7

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    %---- classic real 6dBperoctave peakphased and delayed ----    
    % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
    [b,a] = gammatone(fc,fs,'classic','real','peakphase_new','6dBperoctave');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);
    
    % Find peak at envelope maximum
    envmax1 = zeros(1,nchannels);
    for ii = 1:nchannels
        % Envelope maximum per channel
        envmax1(ii) = find(abs(outsig1(ii,:)) == max(abs(outsig1(ii,:))));
    end;
        
    % Delay signal
    for ii = 1:nchannels
        % Time to delay
        delay =  zeros(1,max(envmax1)+1 - envmax1(ii));
        % Add delay
        outsig1(ii,:) = [delay outsig1(ii,1:N - length(delay))]; 
    end;
    
    % Type of implemantion for headline;
    type1   = 'classic / peakphased(new) / real / +/-6dB / delayed'; 

    % Plot 
    figure ('units','normalized','outerposition',[0 0.05 1 0.475])
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,8*real(outsig1(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, -4, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
    %---- allpole complex 6dBperoctave peakphased and delayed ----
    % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
    [b,a] = gammatone(fc,fs,'allpole','complex','peakphase_new','6dBperoctave');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Find peak at envelope maximum
    envmax2 = zeros(1,nchannels);
    for ii = 1:nchannels
        % Envelope maximum per channel
        envmax2(ii) = find(abs(outsig2(ii,:)) == max(abs(outsig2(ii,:))));
    end;
        
    % Delay signal
    for ii = 1:nchannels
        % Time to delay
        delay =  zeros(1,max(envmax2)+1 - envmax2(ii));
        % Add delay
        outsig2(ii,:) = [delay outsig2(ii,1:N - length(delay))]; 
    end;
    
    % Type of implemantion for headline;
    type2   = 'allpole / peakphased(new) / complex / +/-6dB / delayed';   
    
    % Plot;    
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,8*real(outsig2(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type2)])        
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, -4, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    