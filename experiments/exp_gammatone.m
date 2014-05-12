function exp_gammatone(varargin)
%   Usage: exp_gammatone(flags)
%
%   `exp_gammatone(flags)` reproduces some figures from the papers
%
%  "An efficient auditory filterbank based on the gammatone function"
%  by Roy Patterson and Ian Nimmo-Smith and John Holdsworth and Peter Rice,
%  Annex B of the svos final report,1987.
%
%  "All-Pole models of auditory filtering" by R.F.Lyon,
%  published in "Diversity in Auditory Mechanics", Lewis et al.(eds),
%  World Scientific Publishing, Singapore, 1997, pp. 205-211.
%
%  "Frequency analysis and synthesis using a gammatone filterbank"
%  by V. Hohmann, Acta Acustica united with Acustica Vol.88, 2002
%
%
%   The following flags can be specified;
%
%     'fig3patterson1987'
%
%                   Reproduce Fig.3 similiar to Patterson,1987 efficient:
%                   A cochleagram of four cycles of the [ae] in "past"
%                   produced by by a gammatone filterbank without phase 
%                   compensation. The triangular objects are the upper
%                   three formants of the vowel. The duration of each
%                   period is 8 ms....!!!!
%                       
%
%     'fig5patterson1987'
%
%                   Reproduce Fig.5 from Patterson,1987 efficient:
%                   An array of gamma impulse responses for a 24-channel
%                   auditory filterbank.(lower portion), and the equvalent
%                   array of gamma envelopes (upper portion). The range of
%                   abcissa is 25 ms; the filter centre frequencies range 
%                   from 100 - 4,000 Hz.
%
%                   
%     'fig10apatterson1987'
%                   
%                   Reproduce Fig.10a from Patterson, 1987 eficcient:
%                   The impulse responses for a gammatone auditory filterbank
%                   without phase compensation. The filterbank contains 37 
%                   channels ranging from 100 to 5,000 Hz. The range of abcissa
%                   is 25 ms.
%
%
%     'fig10cpatterson1987'
%                   
%                   Reproduce Fig.10c from Patterson, 1987 efficient:
%                   The impulse responses for a gammatone auditory filterbank
%                   with envelope and fine structure phase-compensation; that is,
%                   the envelope peaks have been aligned and then a fine structure
%                   peak has been aligned with the envelope peak. The filterbank
%                   contains 37 channels ranging from 100 to 5,000 Hz.
%                   The range of abcissa is 25 ms.
%
%
%     'fig11patterson1987'
%
%                   Reproduce similiar to Fig.11 from Patterson, 1987 efficient:
%                   The impulse responses of a gammatone filterbank from a
%                   pulsetrain.
%
%
%     'fig1lyon1997'
%
%                   Reproduce parts of Fig.1 from Lyon, 1997 (no DAPGF):
%                   Comparison of GTF and APGF transfer function for two different
%                   values of the real part of the pole position. Notice that
%                   the ordering of gain near the peak for the classic gammatone
%                   filter is not maintained in the tail. The variation in the
%                   GTF tail is not accounted for by the usual phase-independent
%                   symetric GTF approximation.
%
%
%     'fig2lyon1997'
%
%                   Reproduce parts of Fig.2 from Lyon, 1997 (no DAPGF):
%                   Impulse responses of the APGF and sine-phased GTF from 
%                   Figure 1. Note that the GTF's zero crossings are equally
%                   spaced in time, while those of the APGF are stretched out
%                   in early cycles.
%
%
%     'fig1hohmann2002'
%
%                   Reproduce Fig.l from Hohmann, 2002:
%                   Impulse response of the example Gammatone filter
%                   (center frequency fc = 1000 Hz; 3-db bandwidth fb = 100 Hz;
%                   sampling frequency fs = 10kHz). Solid and dashed lines
%                   show the real and imaginary part of the filter output,
%                   respectively. The absolute value of the filter output
%                   (dashdotted line) clearly represents the envelope.
%
%
%     'fig2hohmann2002'
%
%                   Reproduce Fig.2 from Hohmann, 2002:
%                   Frequency response of the example Gammatone filter
%                   (upper two panels) and of the real-to-imaginary
%                   response(lower two panels). Pi/2 was added to the phase
%                   of the latter (see text). The frequency axis goes up to
%                   half the sampling rate (z=pi).
%       
%
%     'fig3hohmann2002'
%
%                   Reproduce Fig.3 from Hohmann, 2002:
%                   Magnitude frequency response of the Gammatone
%                   filterbank. In this example, the filter channel density
%                   is 1 on the ERB scale and the filter bandwidth is 1
%                   ERBaud. The sampling frequency was 16276Hz and the
%                   lower and upper boundary for the center frequencies
%                   were 70Hz and 6.7kHz, respectively.
%
%
%     'fig4hohmann2002'
%
%                   Treatment of an impulse response with envelope maximum
%                   to the left of the desired group delay (at sample 65).
%                   The original complex impulse response (real part and
%                   envelope plotted in the upper panel) is multiplied with
%                   a complex factor and delayed so that the envelope maximum 
%                   and the maximum of the real part coincide with the desired
%                   group delay (lower panel).
%
%
%   Examples:
%   ---------
%
%  To display Fig. 3 Patterson, 1987 efficient use :::
%
%     exp_gammatone('fig3patterson1987');
%
%   To display Fig. 5 Patterson, 1987 efficient use :::
%
%     exp_gammatone('fig5patterson1987');
%
%  To display Fig. 10a Patterson, 1987 efficient at use :::
%
%     exp_gammatone('fig10apatterson1987');
%
%  To display Fig. 11 Patterson, 1987 efficient use :::
%
%     exp_gammatone('fig11patterson1987');
%
%  To display Fig. 1 Lyon, 1997 use :::
%
%     exp_gammatone('fig1lyon1997');
%
%  To display Fig. 2 Lyon, 1997 use :::
%
%     exp_gammatone('fig2lyon1997');
%
%  To display Fig. 1 Hohmann, 2002 use :::
%
%     exp_gammatone('fig1hohmann2002');
%
%  To display Fig. 2 Hohmann, 2002 use :::
%
%     exp_gammatone('fig2hohmann2002');
%
%  To display Fig. 3 Hohmann, 2002 use :::
%
%     exp_gammatone('fig3hohmann2002');
%
%  To display Fig. 4 Hohmann, 2002 use :::
%
%     exp_gammatone('fig4hohmann2002');
%
%
% AUTHOR: Christian Klemenschitz, 2014

%% ------ Check input options --------------------------------------------
    definput.import={'amtredofile'};
    definput.keyvals.FontSize = 12;
    definput.keyvals.MarkerSize = 6;
    definput.flags.type = {'missingflag', 'fig3patterson1987'...
    'fig5patterson1987' ,'fig10apatterson1987', 'fig10cpatterson1987', ...
    'fig11patterson1987', 'fig1lyon1997', 'fig2lyon1997',...
    'fig1hohmann2002', 'fig2hohmann2002', 'fig3hohmann2002', 'fig4hohmann2002'};

    % Parse input options
    [flags,~]  = ltfatarghelper({'FontSize','MarkerSize'},definput,varargin);

    if flags.do_missingflag
        flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
                   sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
        error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
    end;

%% Pattersons Paper - Classic Gammatone Filter
% Figure 3 
% gammatone 'classic' (casualphase, real)

    if flags.do_fig3patterson1987
        % Input parameters
        %load('a_from_past_us.mat'); % Load mat-file with vowel a in variable insig and sampling frequency in variable fs;
        load('a_from_past_uk.mat'); % Load mat-file with vowel a in variable insig and sampling frequency in variable fs;
        ts = 1/fs;                  % Time between sampling points in s;
        insig = insig(:,1).';       % Extract one channel from signal;
        nchannels = 189;            % Number of channels in filterbank
        N=length(insig);            % Length of signal;
        treal = (1:N)*ts*1000;      % Time axis in ms
        flow = 100;                 % Lowest center frequency;
        fhigh = 4000;               % Highest center frequency;
        fc = zeros(1,nchannels);    % Vector with 37 equally spaced center frequencies;
        bw = ((fhigh-flow)/nchannels); % Distance between frequency channels;
        fc(1) = flow + bw;
        for ii = 2:189
            fc(ii) = fc(ii-1) + bw;
        end;
        
        % Derives filter coefficients for erb spaced channels
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Plot
        figure 
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(treal,12*outsig(ii,:)+dy)
           dy = dy+1;
        end;
        title '189 channel array of gammatone impulse responses of vowel a from past'
        xlabel 'Time (ms)'
        ylabel '189 equal spaced channels by 20.6 Hz from ~ 120.6Hz - 4000Hz'
        set(gca, 'Xlim', [19 80], 'YLim',[0 189])
        hold off
    end;
    
%% Pattersons Paper - Classic Gammatone Filter
% Figure 5 
% gammatone 'classic' (casualphase, real)

    if flags.do_fig5patterson1987
        % Input parameters
        fs = 8000;              % Sampling frequency;
        ts = 1/fs;              % Time between sampling points in s
        N=4096;                 % Length of signal
        treal = (1:N)*ts*1000;  % Time axis in ms
        flow = 100;             % ERB lowest center frequency;
        fhigh = 4000;           % ERB highest center frequency;
        insig = (1:N);          % Impulse as input signal
        insig(:) = 0;
        insig(1) = 1;
        
        % Derives filter coefficients for erb spaced channels
        [b,a] = gammatone(erbspacebw(flow,fhigh),fs,'classic');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Lower figure
        figure('units','normalized','outerposition',[0.5 0.1 0.5 0.9])
        subplot(2,1,2)
        hold on
        dy = 0;
        for ii = 1:size(outsig,1)
           plot(treal,5*outsig(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title '24 channel array of gamma impulse responses (classic)'
        xlabel 'Time (ms)'
        ylabel '24 erb-spaced channels from 100Hz - 4000Hz'
        set(gca, 'XLim',[0 25],'YLim',[0 24])
        hold off
        
        % Upper figure
        subplot(2,1,1)
        hold on
        dy = 0;
        for ii = 1:size(outsig,1)
           env = abs(hilbert(outsig(ii,:),N));
           plot(treal,(5*env)+dy)
           dy = dy+1;
        end;
        title '24 channel array of gamma impulse responses envelopes (classic)'
        xlabel 'Time (ms)'
        ylabel '24 erb-spaced channels from 100Hz - 4000Hz'
        set(gca, 'XLim',[0 25],'YLim',[0 24.5])
        hold off
    end;
    
%% Pattersons Paper 
% Figure 10a
% gammatone 'classic' (casualphase, real)

    if flags.do_fig10apatterson1987
        % Input parameters
        fs = 10000;                 % Sampling frequency;
        ts = 1/fs;                  % Time between sampling points in s;
        N=4096;                     % Length of signal;
        treal = (1:N)*ts*1000;      % Time axis in ms;
        nchannels = 37;             % Number of channels in filterbank;
        flow = 100;                 % ERB lowest center frequency in Hz;
        fhigh = 4000;               % ERB highest center frequency in Hz;
        fc = zeros(1,nchannels);    % Vector with 37 equally spaced center frequencies;
        bw = ((fhigh-flow)/nchannels); % Distance between frequency channels;
        fc(1) = bw;
        for i = 2:nchannels
            fc(ii) = fc(ii-1) + bw;
        end;
        insig = zeros(1,N);         % Impulse as input signal 
        insig(1) = 1;
        
        % Derives filter coefficients for 37 equally spaced channels
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Plot
        figure
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(treal,4*outsig(ii,:) + dy)
           dy = dy+1;
        end;
        clear dy;
        title 'Gammatone impulse response filterbank without phase-compensation'
        xlabel 'Time (ms)'
        ylabel '37 equal spaced channels by 132,43Hz from ~ 135Hz - 5000Hz'
        set(gca, 'XLim',[0 25],'YLim',[0 37])
        hold off
    end;

%% Pattersons Paper
% Figure 10c
% gammatone 'classic','peakphase','complex'

    if flags.do_fig10cpat
        % Input parameters
        fs = 10000;                 % Sampling frequency;
        ts = 1/fs;                  % Time between sampling points in s
        N=4096;                     % Length of signal
        treal = (1:N)*ts*1000;      % Time axis in ms
        nchannels = 37;             % Number of channels in filterbank;
        flow = 100;                 % ERB lowest center frequency in Hz;
        fhigh = 5000;               % Highest center frequency
        fc = zeros(1,nchannels);    % Vector with 37 equally spaced center frequencies;
        bw = (fhigh-flow)/nchannels; % Distance between frequency channels;
        fc(1) = bw;
        for ii = 2:nchannels
            fc(ii) = fc(ii-1) + bw;
        end;
        insig = zeros(1,N);     % Impulse as input signal 
        insig(1) = 1;

        % Derives filter coefficients for 37 equally spaced channels
        [b,a] = gammatone(fc,fs,'classic', 'peakphase');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Find peak at envelope maximum
        outenv = zeros(37,N);      
        envmax = zeros(1,37);
        for ii = 1:nchannels
            % Envelope of filtered impulse signal per channel
            outenv(ii,:) = hilbert(real(outsig(ii,:)),N);
            % Envelope maximum per channel
            envmax(ii) = find(abs(outenv(ii,:)) == max(abs(outenv(ii,:))));
        end;
        
        % Adding Delay
        desiredpeak  = max(envmax)+1;
        outsigdelay = zeros(37,N);
        for ii = 1:nchannels
            % Time to delay
            delay =  zeros(1,desiredpeak - envmax(ii));
            % Add delay
            outsigdelay(ii,:) = [delay outsig(ii,1:size(outsig,2) - length(delay))];
        end;
        
        %Plot
        figure
        hold on
        dy = 0;
        for ii = 1:nchannels
           plot(treal,8*real(outsigdelay(ii,:))+dy)
           dy = dy+1;
        end;
        clear dy;
        title 'Gammatone filterbank with envelope phase-compensation and delayed channel'
        xlabel 'Time (ms)'
        ylabel '37 equal spaced channels by 132,43Hz from ~ 135Hz - 5000Hz'
        set(gca, 'XLim', [0 25], 'YLim',[0 37])
        hold off
    end;
    
%% Pattersons Paper 
% Figure 11 
% gammatone 'classic' (casualphase, real) pulsetrain
  
    if flags.do_fig11patterson1987
        % Input parameters
        fs = 10000;                 % Sampling frequency;
        ts = 1/fs;                  % Time between sampling points in s;
        N=4096;                     % Length of signal;
        treal = (1:N)*ts*1000;      % Time axis in ms;
        nchannels = 37;             % Number of channels in filterbank;
        flow = 100;                 % ERB lowest center frequency in Hz;
        fhigh = 4000;               % ERB highest center frequency in Hz;
        fc = zeros(1,nchannels);    % Vector with 37 equally spaced center frequencies;
        bw = ((fhigh-flow)/nchannels); % Distance between frequency channels;
        fc(1) = bw;
        for ii = 2:nchannels
            fc(ii) = fc(ii-1) + bw;
        end;
        insig = zeros(1,N);     % Impulse as input signal 
        insig(160) = 1;         % at 16 ms,
        insig(240) = 1;         % at 24 ms,
        insig(320) = 1;         % at 32 ms;

        % Derives filter coefficients for erb spaced channels
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
              
        % Plot
        figure('units','normalized','outerposition',[0.5 0.1 0.5 0.9])
        hta1 = subplot(3,1,1);
        posAxes1 = get(hta1, 'Position');
        posAxes1(2) = posAxes1(2) + (1 - posAxes1(2))*5/8;
        posAxes1(4) = posAxes1(4)/4;
        set(hta1, 'Position',posAxes1); 
        plot(treal,insig);
        title 'Pulsetrain'
        ylabel 'Amplitude'
        set(gca, 'XLim',[14.5 39],'YLim',[0 1])
        dy = 0;
        hta2 = subplot(3,1,2);
        posAxes2 = get(hta2, 'Position');
        posAxes2(2) = posAxes1(2)- posAxes2(2);
        posAxes2(4) = posAxes1(4)*6;
        set(hta2, 'Position',posAxes2);
        hold on
        for ii = 1:size(outsig,1)
           plot(treal,4*outsig(ii,:)+dy)
           dy = dy+1;
        end;
        title 'Gammatone impulse response filterbank but without envelope phase-compensation'
        ylabel '37 equal spaced channels'
        set(gca, 'XLim',[14.5 39],'YLim',[0 37])
        hta3 = subplot(3,1,3);
        posAxes3 = get(hta3, 'Position');
        posAxes3(2) =  posAxes3(2) - (1 - posAxes1(2))*1/8; % *3/8
        posAxes3(4) = posAxes1(4)*6;
        set(hta3, 'Position',posAxes3);
 %       title 'test'
        xlabel 'Time (ms)'
        ylabel '37 equal spaced channels'
        set(gca, 'XLim',[12 36],'YLim',[0 37])
        hold off
    end
    
%% Lyons Paper - Allpass Gammatone Filter
% Figure 1
% gammatone (allpass,casualphase,real')
% gammatone 'classic' (casualphase, real)

    if flags.do_fig1lyon1997
        warning(['FIXME: The real-valued allpass filters produce good results, ' ... 
                 'though there is a fixme for real-valued filters, ' ...
                 'complex-valued filters do not.']);

        flow = 100;
        fhigh = 4000;
        fs = 25000;
        erb = erbspacebw(flow,fhigh);

        n = 6;
        betamul = 3;
        [b,a] = gammatone(erb,fs,n,betamul);
        [h,w] = freqz(b(13,:), a(13,:));
        h = h(:)/h(1);
        nor = real(max(h));

        figure
        semilogx(w/(erb(13)/fs*2*pi), 20*log10(abs(h)),'r-');
        title 'Transfer functions of classic and allpass gammatone filters for b=omeg_r/3 and b=omeg_r/2'
        xlabel 'omeg/omeg_r'
        ylabel 'Magnitude gain (db)'
        set(gca,'XLim',[0.03 3],'YLim',[-40, 40]) 
        hold on

        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        [h,w] = freqz(b(13,:), a(13,:));
        h=abs(h) * nor;
        semilogx(w/(erb(13)/fs*2*pi), 20*log10(abs(h)),'b--');

        betamul = 2;
        [b,a] = gammatone(erb,fs,n,betamul);
        [h,w] = freqz(b(13,:), a(13,:));
        h = h(:)/h(1);
        nor = real(max(h));
        semilogx(w/(erb(13)/fs*2*pi), 20*log10(abs(h)),'r-');

        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        [h,w] = freqz(b(13,:), a(13,:));
        h=abs(h) * nor;
        semilogx(w/(erb(13)/fs*2*pi), 20*log10(abs(h)),'b--');
        hold off
    end;
%% Lyons Paper 
% Figure 2
% gammatone 'classic' (casualphase, real)
% gammatone (allpass,casualphase) 'complex'

    if flags.do_fig2lyon1997

        flow =100;
        fhigh = 4000;
        fs = 25000;
        erb = erbspacebw(flow,fhigh);
        N = 2048;
        insig = (1:N);
        insig(:) = 0;
        insig(1024) = 1;

        n = 6;
        betamul = 3;
        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
   
        figure
        subplot(2,1,1)
        plot(outsig(13,:),'b--')
        title 'Impulse responses of classic and allpass gammatone filters for b=omeg_r/3'
        xlabel 'Sample points'
        ylabel 'Amplitude'
        set(gca, 'XLim',[1020 1084],'YLim',[-0.25 0.25])
        hold on
        
        %[b,a] = gammatone(erb,fs,n,betamul);
        [b,a] = gammatone(erb,fs,n,betamul,'complex');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        plot(outsig(13,:),'r-')
        hold off

        n = 6;
        betamul = 2;
        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);

        subplot(2,1,2)
        plot(outsig(13,:),'b--')
        title 'Impulse responses of classic and allpass gammatone filters for b=omeg_r/2'
        xlabel 'Sample points'
        ylabel 'Amplitude'
        set(gca, 'XLim',[1020 1084],'YLim',[-0.25 0.25])
        hold on

        %[b,a] = gammatone(erb,fs,n,betamul);
        [b,a] = gammatone(erb,fs,n,betamul,'complex');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        plot(outsig(13,:),'r-')
        hold off
    end;
%% Hohmanns paper
% Figure 1
% gammatone (allpass,casualphase) 'complex'

    if flags.do_fig1hohmann2002
        fs = 25000;
        fc = 1000;
        N = 4096;
        insig = (1:N);
        insig(:) = 0;
        insig(1) = 1;

        [b,a] = gammatone(fc,fs,'complex');
        outsig = (ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        outsigenv = hilbert(outsig,fs);

        figure  
        plot(abs(outsigenv),'b-.')
        set(gca, 'XLim',[0 200],'YLim',[-0.04 0.04])
        hold on
        plot(real(outsig),'r-')
        plot(imag(outsig),'g--')
        hold off
    end;
    
%% Hohmanns paper
% Figure 2
% gammatone (allpass,casualphase) 'complex'

    if flags.do_fig2hohmann2002
        fc = 1000;
        fs = 25000;
        fn = 5000;
        Tn = 1/fn;
        N = 4096;
        df = fs/N;
        xfn = (0:df:fn-df)*Tn;
        insig = (1:N);
        insig(:) = 0;
        insig(1) = 1;

        [b,a] = gammatone(fc,fs,'complex');
        outsig = (ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Upper two panels
        outfftr = fft(2*real(outsig));
        [phi] = phasez(b,a,2048);
        phi = phi+2*pi;
        phi = phi.';
        
        % Lower two panels
        % The real-to-imaginary transfer function was calculated by
        % dividing the FFT-spectra of imaginary and real part of the
        % impulse response.
        outffti = fft(2*imag(outsig));
        outsigrtoi = outffti./outfftr;
        theta = angle(outsigrtoi);
        phirtoi = theta+ pi/2;
               
        figure
        subplot(4,1,1)
        plot(xfn, 20*log10(abs(outfftr(1:N/2))))
        set(gca,'Xlim',[0 1], 'YLim',[-70 0])
        ylabel('Magnitude/db')
        subplot(4,1,2)
        plot(xfn,phi)
        ylabel('Phase/rad')
        set(gca,'Xlim',[0 1],'Ylim',[-5 5])
        subplot(4,1,3)
        plot(xfn, 20*log10(abs(outsigrtoi(1:N/2))))
        set(gca,'Xlim',[0 1],'Ylim',[-10 10])
        ylabel('Magnitude/db')
        subplot(4,1,4)
        plot(xfn, phirtoi(1:N/2))
        set(gca,'Xlim',[0 1],'Ylim',[-2 2])
        xlabel('Frequency / pi')
        ylabel('Phase/rad')
        grid
    end;
    
%% Hohmanns paper
% Figure 3
% gammatone (allpass,casualphase) 'complex'

     if flags.do_fig3hohmann2002
        flow = 70;
        fhigh = 6700;
        fs = 16276;
        fn = fs/2;
        N = 2048;
        erb = erbspacebw(flow,fhigh);

        [b,a] = gammatone(erb,fs,'complex');
        h = zeros(length(erb),N/4);
        w = zeros(length(erb),N/4);
        
        for ii=1:length(erb)
           [h(ii,:),w(ii,:)] = freqz(b(ii,:),a(ii,:));
        end;
        h(:) = abs(h(:));
        
        figure
        hold on
        for ii = 1:length(erb)
            plot(w(ii,:)/(2*pi)*fs, 20*log10(h(ii,:)))
        end;
        xlabel('Frequency / Hz')
        ylabel('Level / dB')
        set(gca,'Xlim',[0 fn],'Ylim',[-40 0])
        hold off
     end;
     
%% Hohmanns paper
% Figure 4
% gammatone (allpass,casualphase) 'complex'
% Result is delayed and peakphased.
% Case 2a: envelope maximum before desired peak;

    if flags.do_fig4hohmann2002
        flow = 70;
        fhigh = 6700;
        fs = 16276;
        insig = (1:fs);
        insig(:) = 0;
        insig(1) = 1;
        erb = erbspacebw(flow,fhigh);
       
        % The envelope maximum envmax is earlier in time than the desired
        % peak at channels 13 to 30.
        channel = 20;
        desiredpeak = 65;

        [b,a] = gammatone(erb,fs,'complex');
        outsig = (ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);

        outenv = hilbert(real(outsig(channel,:)),fs);
        envmax = find(abs(outenv(1,:)) == max(abs(outenv(1,:))));
        sigmax = find(real(outsig(channel,:)) == max(real(outsig(channel,:))));

        % Equation 18 Phi_k
        phi = erb(channel)*(-2*pi)*(envmax - sigmax)/fs;
       
        % Equation 19
        outsignew = real(outsig(channel,:)) * cos(phi) - imag(outsig(channel,:)) * sin(phi);
          
          % For testing purpose 
%         outenvnew = hilbert(real(outsignew(1,:)),fs);
%         envmaxnew = find(abs(outenvnew(1,:)) == max(abs(outenvnew(1,:))));
%         sigmaxnew = find(real(outsignew(1,:)) == max(real(outsignew(1,:))));
%         
%         warning(['FIXME: The maximum of the real part of the impulse response ' ... 
%                  'misses the maximum of the envelope in channel '  num2str(channel)...
%                  ' by ' num2str(envmaxnew - sigmaxnew) ' sample(s).']);
                
        % Equation 20
        delay = outsignew(end-desiredpeak+envmax+1:end);
        
        % Equation 21
        outsigdelay= [delay outsignew(1:length(outsignew) - delay)];
        outsigdelay = outsigdelay(1:length(outsignew) - delay);
        
        outsigdelayenv = hilbert(real(outsigdelay),fs);
        
        figure
        subplot(2,1,1)
        plot(real(outsig(channel,:)),'b-')
        hold on
        set(gca, 'XLim',[0 200],'YLim',[-0.04 0.04])
        xlabel('Sample')
        ylabel('Amplitude')
        plot(abs(outenv),'r-')
        hold off
        subplot(2,1,2)
        plot (real(outsigdelay),'b-')
        hold on
        set(gca, 'XLim',[0 200],'YLim',[-0.04 0.04])
        xlabel('Sample')
        ylabel('Amplitude')
        plot(abs(outsigdelayenv),'r-')
        hold off
    end;
     
end