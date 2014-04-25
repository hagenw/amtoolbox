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
%     'fig5pat'     Reproduce Fig.5 from Patterson,1987 efficient:
%                   An array of gamma impulse responses for a 24-channel
%                   auditory filterbank.(lower portion), and the equvalent
%                   array of gamma envelopes (upper portion). The range of
%                   abcissa is 25 ms; the filter centre frequencies range 
%                   from 100 - 4,000 Hz.
%                   
%     'fig10apat'   Reproduce Fig.10a from Patterson, 1987 eficcient:
%                   The impulse responses for a gammatone auditory filterbank
%                   without phase compensation. The filterbank contains 37 
%                   channels ranging from 100 to 5,000 Hz. The range of abcissa
%                   is 25 ms.
%               
%     'fig10cpat'   Reproduce Fig.10c from Patterson, 1987 efficient:
%                   The impulse responses for a gammatone auditory filterbank
%                   with envelope and fine structure phase-compensation; that is,
%                   the envelope peaks have been aligned and then a fine structure
%                   peak has been aligned with the envelope peak. The filterbank
%                   contains 37 channels ranging from 100 to 5,000 Hz.
%                   The range of abcissa is 25 ms.
%
%     'fig11pat'    Reproduce similiar to Fig.11 from Patterson, 1987 efficient:
%                   The impulse responses of a gammatone filterbank from a
%                   pulsetrain.
%
%     'fig1ly'      Reproduce parts of Fig.1 from Lyon, 1997 (no DAPGF):
%                   Comparison of GTF and APGF transfer function for two different
%                   values of the real part of the pole position. Notice that
%                   the ordering of gain near the peak for the classic gammatone
%                   filter is not maintained in the tail. The variation in the
%                   GTF tail is not accounted for by the usual phase-independent
%                   symetric GTF approximation.
%
%     'fig2ly'      Reproduce parts of Fig.2 from Lyon, 1997 (no DAPGF):
%                   Impulse responses of the APGF and sine-phased GTF from 
%                   Figure 1. Note that the GTF's zero crossings are equally
%                   spaced in time, while those of the APGF are stretched out
%                   in early cycles.
%               
%     'fig1ho'      Reproduce parts of Fig.l from Hohmann, 2002:
%                   Impulse response of the example Gammatone filter
%                   (center frequency fc = 1000 Hz; 3-db bandwidth fb = 100 Hz;
%                   sampling frequency fs = 10kHz). The solid and dashed lines
%                   show the real and imaginary part of the filter output,
%                   respectively. The absolute value of the filter output
%                   (dashdotted line) clearly represents the envelope.
%
%     'fig2ho'      Reproduce Fig.2 from Hohmann, 2002:
%                   Still in Progress
%
%   Examples:
%   ---------
%
%   To display Fig. 5 Patterson, 1987 efficient use :::
%
%     exp_gammatone('fig5pat');
%
%  To display Fig. 10a Patterson, 1987 efficient at use :::
%
%     exp_gammatone('fig10apat');
%
%  To display Fig. 10c Patterson, 1987 efficient use :::
%
%     exp_gammatone('fig10cpat');
%
%  To display Fig. 11 Patterson, 1987 efficient use :::
%
%     exp_gammatone('fig11pat');
%
%  To display Fig. 1 Lyon, 1997 use :::
%
%     exp_gammatone('fig1ly');
%
%  To display Fig. 2 Lyon, 1997 use :::
%
%     exp_gammatone('fig2ly');
%
%  To display Fig. 1 Hohmann, 2002 use :::
%
%     exp_gammatone('fig1ho');
%
%  To display Fig. 2 Hohmann, 2002 use :::
%
%     exp_gammatone('fig2ho');
%
%
% AUTHOR: Christian Klemenschitz, 2014

%% ------ Check input options --------------------------------------------
    definput.import={'amtredofile'};
    definput.keyvals.FontSize = 12;
    definput.keyvals.MarkerSize = 6;
    definput.flags.type = {'missingflag',...
    'fig5pat' ,'fig10apat', 'fig10cpat', 'fig11pat', 'fig1ly', 'fig2ly',...
    'fig1ho', 'fig2ho'};

    % Parse input options
    [flags,~]  = ltfatarghelper({'FontSize','MarkerSize'},definput,varargin);

    if flags.do_missingflag
        flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
                   sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
        error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
    end;

%% Pattersons Paper - Classic Gammatone Filter
% Figure 5 
% gammatone 'classic' (casualphase, real)

    if flags.do_fig5pat
        flow=100;
        fhigh=4000;
        fs=8000;
        N=2048;
        insig=(1:N);
        insig(:)=0;
        insig(1)=1;

        [b,a] = gammatone(erbspacebw(flow,fhigh),fs,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);

        figure('units','normalized','outerposition',[0.5 0.1 0.5 0.9])
        subplot(2,1,2)
        hold on
        dy = 0;
        for i=1:size(outsig,1)
           plot(outsig(i,:)+dy)
           dy = dy+0.2;
        end
        clear dy;
        title '24 channel array of gamma impulse responses (classic)'
        xlabel 'Sample points of 25ms (fs=8000Hz)'
        ylabel '24 erb-spaced channels from 100Hz - 4000Hz'
        set(gca, 'XLim',[0 200],'YLim',[0 4.8])
        hold off

        subplot(2,1,1)
        hold on
        dy = 0;
        for ii=1:size(outsig,1)
           env= abs(hilbert(outsig(ii,:),fs));
           plot(env+dy)
           dy = dy+0.2;
        end
        title '24 channel array of gamma impulse responses envelopes (classic)'
        xlabel 'Sample points of 25ms (fs=8000Hz)'
        ylabel '24 erb-spaced channels from 100Hz - 4000Hz'
        set(gca, 'XLim',[0 200],'YLim',[0 5])
        hold off


    end
%% Figure 10a
% gammatone 'classic' (casualphase, real)

    if flags.do_fig10apat
        fhigh=5000;
        fs=10000;
        bw=zeros(1,37);
        bw(1)=(fhigh/37)-17.57;
        N=2048;
        insig=(1:N);
        insig(:)=0;
        insig(1024)=1;
        for i=2:37
        bw(i)=bw(i-1)+(fhigh/37);
        end

        [b,a] = gammatone(bw,fs,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);

        figure
        hold on
        dy = 0;
        for i=1:size(outsig,1)
           plot((outsig(i,:))+dy)
           dy = dy+0.1;
        end
        clear dy;
        title 'Gammatone impulse response filterbank without phase-compensation'
        xlabel 'Sample points of ~14.3ms (fs=10000Hz)'
        ylabel '37 equal spaced channels by 132,43Hz from ~ 100Hz - 5000Hz'
        set(gca, 'XLim',[1022 1094],'YLim',[0 3.7])
        hold off
    end
%% Figure 10c
% gammatone 'classic','peakphase','complex'

    if flags.do_fig10cpat

        fhigh=5000;
        fs=10000;
        bw=zeros(1,37);
        bw(1)=(fhigh/37)-17.57;
        N=2048;
        insig=(1:N);
        insig(:)=0;
        insig(1024)=1;
        for i=2:37
        bw(i)=bw(i-1)+(fhigh/37);
        end

        [b,a] = gammatone(bw,fs,'classic','peakphase','complex');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);

        figure
        hold on
        dy = 0;
        for i=1:size(outsig,1)
           plot((outsig(i,:))+dy)
           dy = dy+0.1;
        end
        clear dy;
        title 'Gammatone filterbank with envelope and fine-structure phase-compensation'
        xlabel 'Sample points of ~14.3ms (fs=10000Hz)'
        ylabel '37 equal spaced channels by 132,43 Hz // ~ 100 Hz - 5000 Hz'
        set(gca, 'XLim',[1022 1094],'YLim',[0 3.7])
        hold off
    end  
%% Figure 11 
% gammatone 'classic' (casualphase, real) pulsetrain
  
    if flags.do_fig11pat
        fhigh=5000;
        fs=10000;
        bw=zeros(1,37);
        bw(1)=(fhigh/37)-17.57;
        N=2048;
        insig=(1:N);
        insig(:)=0;
        insig(1)=1;
        insig(40)=1;
        insig(80)=1;

        for i=2:37
        bw(i)=bw(i-1)+(fhigh/37);
        end

        [b,a] = gammatone(bw,fs,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);

        figure('units','normalized','outerposition',[0.5 0.1 0.5 0.9])
        subplot(2,1,1)
        plot(insig)
        title 'Pulsetrain'
        xlabel 'Sample points of 12ms (fs=10000Hz)'
        ylabel 'Amplitude'
        set(gca, 'XLim',[0 120],'YLim',[0 1])
        dy = 0;
        subplot(2,1,2)
        hold on
        for i=1:size(outsig,1)
           plot((outsig(i,:))+dy)
           dy = dy+0.1;
        end
        clear dy;
        title 'Gammatone impulse response filterbank without phase-compensation'
        xlabel 'Sample points of 12ms (fs=10000Hz)'
        ylabel '37 equal spaced channels from ~ 100Hz - 5000Hz'
        set(gca, 'XLim',[0 120],'YLim',[0 3.7])
        hold off

    end
%% Lyons Paper - Allpass Gammatone Filter
% Figure 1
% gammatone (allpass,casualphase,real')
% gammatone 'classic' (casualphase, real)

    if flags.do_fig1ly
        warning(['FIXME: The real-valued allpass filters produce good results, ' ... 
                 'though there is a fixme for real-valued filters, ' ...
                 'complex-valued filters do not.']);

        flow=100;
        fhigh=4000;
        fs=8000;
        erb=erbspacebw(flow,fhigh);

        n=6;
        betamul=3;
        [b,a] = gammatone(erb,fs,n,betamul);
        [h,w] = freqz(b(13,:),a(13,:));
        h=h(:)/h(1);
        nor=real(max(h));

        figure
        semilogx(w/(erb(13)/fs*2*pi),20*log10(abs(h)),'r-');
        title 'Transfer functions of classic and allpass gammatone filters for b=omeg_r/3 and b=omeg_r/2'
        xlabel 'omeg/omeg_r'
        ylabel 'Magnitude gain (db)'
        set(gca,'XLim',[0.03 3],'YLim',[-40, 40]) 
        hold on

        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        [h,w] = freqz(b(13,:),a(13,:));
        h=abs(h)*nor;
        semilogx(w/(erb(13)/fs*2*pi),20*log10(abs(h)),'b--');

        betamul=2;
        [b,a] = gammatone(erb,fs,n,betamul);
        [h,w] = freqz(b(13,:),a(13,:));
        h=h(:)/h(1);
        nor=real(max(h));
        semilogx(w/(erb(13)/fs*2*pi),20*log10(abs(h)),'r-');

        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        [h,w] = freqz(b(13,:),a(13,:));
        h=abs(h)*nor;
        semilogx(w/(erb(13)/fs*2*pi),20*log10(abs(h)),'b--');
        hold off
    end
%% Figure 2
% gammatone 'classic' (casualphase, real)
% gammatone (allpass,casualphase) 'complex'

    if flags.do_fig2ly
        warning(['FIXME: The complex-valued allpass filters are not scaled ' ...
                 'correctly too.']);
        flow =100;
        fhigh=4000;
        fs=8000;
        erb=erbspacebw(flow,fhigh);
        N=2048;
        insig=(1:N);
        insig(:)=0;
        insig(1024)=1;

        n=6;
        betamul=3;
        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);
   
        figure
        subplot(2,1,1)
        plot(outsig(13,:),'b--')
        title 'Impulse responses of classic and allpass gammatone filters for b=omeg_r/3'
        xlabel 'Sample points'
        ylabel 'Amplitude'
        set(gca, 'XLim',[1020 1084],'YLim',[-0.25 0.25])
        hold on

        [b,a] = gammatone(erb,fs,n,betamul,'complex');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);
        plot(outsig(13,:),'r-')
        hold off

        n=6;
        betamul=2;
        [b,a] = gammatone(erb,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);

        subplot(2,1,2)
        plot(outsig(13,:),'b--')
        title 'Impulse responses of classic and allpass gammatone filters for b=omeg_r/2'
        xlabel 'Sample points'
        ylabel 'Amplitude'
        set(gca, 'XLim',[1020 1084],'YLim',[-0.25 0.25])
        hold on

        [b,a] = gammatone(erb,fs,n,betamul,'complex');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);
        plot(outsig(13,:),'r-')
        hold off
    end
%% Hohmanns paper
% Figure 1
% gammatone (allpass,casualphase) 'complex'

    if flags.do_fig1ho
%         warning(['FIXME: The complex-valued allpass filters produce ' ...
%                  'no imaginary part.']);
        fs=10000;
        fc=1000;
        N=4096;
        insig=(1:N);
        insig(:)=0;
        insig(1)=1;

        [b,a] = gammatone(fc,fs,'complex');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);
        outsig = hilbert(outsig,fs);

        figure  
        plot(abs(outsig),'b-.')
        set(gca, 'XLim',[0 200],'YLim',[-0.04 0.04])
        hold on
        plot(real(outsig),'r-')
        plot(imag(outsig),'g--')
        hold off
    end
%% Figure 2
%
    if flags.do_fig2ho
         warning(['FIXME: Still in progress..' ...
                 '.']);
        fs=10000;
        fc=1000;
        N=4096;
        insig=(1:N);
        insig(:)=0;
        insig(1)=1;

        [b,a] = gammatone(fc,fs,'complex');
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig=permute(outsig,[3 2 1]);
        outfft= fft(outsig);
        figure
        stem( outfft(N:end), 'b.-')
        %axis([-Fn Fn 0 max_y])
        title('Amplitudengang')
        ylabel('Amplitude')
        xlabel(['Auflösung: ',num2str(N),' Hz Frequenz in Hz'])
        grid
    end
end