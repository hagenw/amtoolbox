function [dataOut] = exp_moore1997(varargin)
%EXP_MOORE2002 Figures several papers by Moore et al.
%

%
%   Usage:
%     [dataOut] = exp_kelvasa2015('fig8a')
%
%   The following flags can be specified;
%     'fig2'    Reproduce Fig. 2 of moore1997.
%
%     'fig3'    Reproduce Fig. 3 of moore1997.
%
%     'fig5'    Reproduce Fig. 5 of moore1997.
%
%     'fig6'    Reproduce Fig. 6 of moore1997.
%
%     'fig7'    Reproduce Fig. 7 of moore1997.
%
%     'fig8'    Reproduce Fig. 8 of moore1997.
%
%


%% Retrieve and compute model paramters
    % Set flags

    definput.flags.type = {'missingflag','fig2','fig3','fig5','fig6','fig7','fig8'};

    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
           flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
           sprintf('%s or %s',definput.flags.type{end-1},...
           definput.flags.type{end})];
           error('%s: You must specify one of the following flags: %s.',...
                 upper(mfilename),flagnames);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig2
        fs = 32000;
        fVec = 20:fs/2;
        data = data_moore2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);
        figure
        semilogx(data.fOuter,data.tfOuter);
        xlim([20,20000])
        ylim([-5,20])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig3
        fs = 32000;
        fVec = 20:fs/2;
        data = data_moore2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);
        figure
        semilogx(data.fMiddle,abs(data.tfMiddle));
        xlim([20,20000])
        ylim([0,40])
        xlabel('Frequency, kHz')
        ylabel('Middle ear effective attenuation')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig5
        % todo: figure not yet the same as in moore1997 paper
        fs = 32000;
        t = linspace(0,1,fs);
        sig = sin(2*pi*1000*t).';
        figure
        for ii=1:9
        inSig = setdbspl(sig,(ii+1)*10);  % plot from 20 to 100 dB in 10 dB steps
        results = moore1997(inSig,fs);
        plot(results.erbN,results.eLdB,'k')
        hold on
        end
        xlim([5,40])
        ylim([0,100])
        xlabel('Number of ERB')
        ylabel('Excitation level, dB')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig6
        erbStep = 0.25;    % according to moore1997, moore2002
        erbFcMin = 50;
        erbFcMax = 15000;
        erbNMin = fc2erbN(erbFcMin);
        erbNMax = fc2erbN(erbFcMax);
        erbN = erbNMin:erbStep:erbNMax;    % numbers of erb bands
        erbFc = erbN2fc(erbN); 
        dataSL = data_moore2002('specLoud','fVec',erbFc);
        gdB = dataSL.g;    % low level gain in cochlea amplifier
        alpha = dataSL.alpha;    % compressive exponent
        figure
        plot(gdB,alpha)
        xlim([-25,0])
        ylim([0.2,0.3])
        xlabel('G, dB')
        ylabel('alpha')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig7
        erbStep = 0.25;    % according to moore1997, moore2002
        erbFcMin = 50;
        erbFcMax = 15000;
        erbNMin = fc2erbN(erbFcMin);
        erbNMax = fc2erbN(erbFcMax);
        erbN = erbNMin:erbStep:erbNMax;    % numbers of erb bands
        erbFc = erbN2fc(erbN); 
        dataSL = data_moore2002('specLoud','fVec',erbFc);
        gdB = dataSL.g;    % low level gain in cochlea amplifier
        a = dataSL.a;    % parameter for linearization around absolute threshold
        figure
        plot(gdB,a)
        xlim([-25,0])
        ylim([4,10])
        xlabel('G, dB')
        ylabel('A')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 8
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig8

    end

end
