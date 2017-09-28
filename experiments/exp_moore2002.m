function [dataOut] = exp_moore2002(varargin)
%EXP_MOORE2002 Figures several papers by Moore et al.
%

%
%   Usage:
%     [dataOut] = exp_kelvasa2015('fig8a')
%
%   The following flags can be specified;
%
%     'fig5'    Reproduce Fig. 5 of moore1997.
%
%     'fig8'    Reproduce Fig. 8 of moore1997.
%
%     'fig1'    Reproduce Fig. 1 of moore2002.
%
%     'fig1b'    Similar to fig1 but using 2006 revised data for middle ear
%                filter.
%
%     'fig2'    Reproduce Fig. 2 of moore2002.
%
%     'fig3'    Reproduce Fig. 3 of moore2002.
%
%     'fig5b'    Reproduce Fig. 5 of moore2002.
%


%% Retrieve and compute model paramters
    % Set flags

    definput.flags.type = {'missingflag','fig5','fig8','fig1','fig1b','fig2',...
                                           'fig3','fig5b'};

    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
           flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
           sprintf('%s or %s',definput.flags.type{end-1},...
           definput.flags.type{end})];
           error('%s: You must specify one of the following flags: %s.',...
                 upper(mfilename),flagnames);
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig5
        fs = 32000;
        t = linspace(0,1,fs);
        inSig = sin(2*pi*1000*t).';
        results = moore2002(inSig,fs);
        semilogx(results.eLdB(500,:))
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 8
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig8

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig1
        fs = 32000;
        fVec = 20:fs/2;
        data = data_moore2002('tfOuterMiddle1997','fieldType','free','fVec',fVec);
        figure
        semilogx(fVec, data.tfOuterMiddle)
        grid on
        xlim([20,16000])
        xlabel('Frequency (Hz)')
        ylabel('Relative Transmission (dB)')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 1b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig1b
        fs = 32000;
        fVec = 20:fs/2;
        data = data_moore2002('tfOuterMiddle2007','fieldType','free','fVec',fVec);
        figure
        semilogx(fVec, data.tfOuterMiddle)
        grid on
        xlim([20,16000])
        xlabel('Frequency (Hz)')
        ylabel('Relative Transmission (dB)')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig2

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig3

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Figure 5b (Fig. 5 of moore2002)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_fig5b

    end

end
