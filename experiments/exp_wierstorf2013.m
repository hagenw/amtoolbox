function output = exp_wierstorf2013(varargin)
%EXP_WIERSTORF2013 Figures from Wierstorf (2013)
%   Usage: output = exp_wierstorf2013(flag)
%
%   `exp_wierstorf2013(flag)` reproduces the results for the figure given
%   by *flag* from the Wierstorf (2013) paper. It will also plot the
%   results.  The format of its output depends on the chosen figure.
%   Note that the original figures in the paper are not plotted with Matlab, but
%   with gnuplot and may look a little bit different.
%   
%   The following flags can be specified;
%
%     'plot'     plot the output of the experiment. This is the default.
%
%     'noplot'   Don't plot, only return data.
%
%     'auto'     Re-calculate the file if it does not exist. Return 1 if the
%                file exist, otherwise 0. This is the default
%
%     'refresh'  Always recalculate the file.
%
%     'cached'   Always use the cached version. Throws an error if the
%                file does not exist.
%
%     'fig1'     Reproduce Fig.1 from Wierstorf (2013). The localization error
%                for a typical stereophony setup is calculated and shown for the
%                whole listening are, sampled with 21x21 point.
%
%     'fig3'     Simulations of the sound field for Wave Field Synthesis for
%                a mono-frequent virtual plane wave with three different
%                frequencies of 1kHz, 2kHz, 5kHz. In addition a spatio-temporal
%                impulse response of the sound field for a broadband plane wave
%                is shown at the time 4.8ms after its start.
%
%     'fig6'     Results from an experiment comparing the localization accuracy
%                for a real point source (loudspeaker) and a simulated point
%                source (binaural synthesis).
%
%     'fig7'     Results from a localization experiment for a virtual point
%                source in Wave Field Synthesis for different positions in the
%                listening area (data are the same as in Fig.10). The line is
%                always starting from a listener position and points towards the
%                direction the listener perceived the auditory event.
%
%     'fig8'     Mapping of ITD values in the first twelve frequency channels to
%                the corresponding azimuth angles in the range -90deg to 90deg.
%                The ITD values are calculated from an HRTF data base with the
%                binaural model after Dietz.
%
%     'fig9'     Deviation of the predicted sound source location with the
%                mapping function from Fig.8 for the same HRTF data set as in
%                Fig.8.
%
%     'fig10'    Results from a localization experiment for a virtual point
%                source in Wave Field Synthesis for different positions in the
%                listening area (data points are the same as in Fig.6). The
%                signals were simulated by binaural synthesis and given also to
%                the Dietz binaural model to predict the localization. The model
%                results are shown as lines. Three different loudspeaker array
%                setups were used.
%
%     'fig11a'   Prediction of the localization for a virtual point source in
%                Wave Field Synthesis in the whole listening area for a linear
%                loudspeaker array.
%
%     'fig11b'   Prediction of the localization for a virtual plane wave in Wave
%                Field Synthesis in the whole listening area for a linear
%                loudspeaker array.
%
%     'fig12a'   Prediction of the localization for a virtual point source in
%                Wave Field Synthesis in the whole listening area for a circular
%                loudspeaker array.
%
%     'fig12b'   Prediction of the localization for a virtual plane wave in Wave
%                Field Synthesis in the whole listening area for a circular
%                loudspeaker array.
%
%   If no flag is given, the function will print the list of valid flags.
%
%   Examples:
%   ---------
%
%   To display Figure 1 use :::
%
%     exp_wierstorf('fig1');
%
%   References: wierstorf2013

%   AUTHOR: Hagen Wierstorf

definput.import={'amtredofile'};
definput.flags.type={'missingflag','fig1','fig3','fig6','fig7','fig8', ...
                    'fig9','fig10','fig11a','fig11b','fig12a','fig12b'};

definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

save_format='-v6';

% Checking for the Sound-Field-Synthesis Toolbox
if ~which('SFS_start')
    error(['%s: you need to install the Sound-Field-Synthesis Toolbox.\n', ...
        'You can download it at https://github.com/sfstoolbox/sfs.\n', ...
        'You need version 0.2.4 of the Toolbox (commit afe5c14359).'], ...
        upper(mfilename));
end
conf = SFS_config;


%% ------ F I G U R E  1 -------------------------------------------------
if flags.do_fig1

    s = [mfilename('fullpath'),'_fig1.mat'];
      
    % listening area
    X = [-2 2];
    Y = [-3.15 -0.15];
    % orientation of the listener (always to the front)
    phi = pi/2;
    % position of the virtual point source
    xs = [0 0];
    src = 'ps';
    % distance between the loudspeaker ion the stereo setup
    L = 2;
  
    if amtredofile(s,flags.redomode)
        [loc_error,aud_event,sound_event,xaxis,yaxis,x0] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'stereo');
        save(save_format,s,'loc_error','aud_event','xaxis','yaxis','x0');
    else
        load(s);
    end;

    output.loc_error = loc_error;
    output.aud_event = aud_event;
    output.xaxis = xaxis;
    output.yaxis = yaxis;
    output.x0 = x0;

    if flags.do_plot
        % ------ Plotting ------
        figure;
        [u,v,~] = pol2cart(rad(aud_event+90),ones(size(aud_event)), ...
            zeros(size(aud_event)));
        quiver(xaxis,yaxis,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0,1);
        xlabel('x/m');
        ylabel('y/m');
    end;

%% ------ F I G U R E  3  ------------------------------------------------
elseif flags.do_fig3

    s = [mfilename('fullpath'),'_fig3.mat'];
 
    % listening area
    X = [-2 2];
    Y = [-2 2];
    % circular array with 56 loudspeakers
    L = 3;
    conf.array = 'circle';
    conf.dx0 = L*pi/56;
    conf.xref = [0 0];
    % plane wave travelling upwards
    xs = [0 1];
    src = 'pw';
    % other neccessary settings
    conf.usebandpass = 0;
    conf.bandpassflow = 0;
    conf.bandpassfhigh = 20000;
    conf.plot.usedb = 0;
    conf.plot.colormap = '';
  
    if amtredofile(s,flags.redomode)
        % get secondary sources and tapering window for plotting
        x0 = secondary_source_positions(L,conf);
        win = [zeros(1,29) ones(1,27)];
        % (a)
        f = 1000;
        [xaxis,yaxis,P_a] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
        % (b)
        f = 2000;
        [~,~,P_b] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
        % (c)
        f = 5000;
        [~,~,P_c] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
        % (d)
        t = 212;
        [~,~,p_d] = wave_field_imp_wfs_25d(X,Y,xs,src,t,L,conf);
        save(save_format,s,'P_a','P_b','P_c','p_d','xaxis','yaxis','x0','win');
    else
        load(s);
    end;

    output.P_a = P_a;
    output.P_b = P_b;
    output.P_c = P_c;
    output.p_d = p_d;
    output.xaxis = xaxis;
    output.yaxis = yaxis;
    output.x0 = x0;
    output.win = win;

    if flags.do_plot
        % ------ Plotting ------
        % (a)
        plot_wavefield(xaxis,yaxis,P_a,x0,win,conf);
        axis([X(1) X(2) Y(1) Y(2)]);
        colorbar;
        xlabel('x/m');
        ylabel('y/m');
        title('(a) f_{pw} = 1kHz');
        % (b)
        plot_wavefield(xaxis,yaxis,P_b,x0,win,conf);
        axis([X(1) X(2) Y(1) Y(2)]);
        colorbar;
        xlabel('x/m');
        ylabel('y/m');
        title('(b) f_{pw} = 2kHz');
        % (c)
        plot_wavefield(xaxis,yaxis,P_c,x0,win,conf);
        axis([X(1) X(2) Y(1) Y(2)]);
        colorbar;
        xlabel('x/m');
        ylabel('y/m');
        title('(c) f_{pw} = 5kHz');
        % (d)
        conf.plot.usedb = 1;
        plot_wavefield(xaxis,yaxis,p_d,x0,win,conf);
        axis([X(1) X(2) Y(1) Y(2)]);
        colorbar;
        xlabel('x/m');
        ylabel('y/m');
        title('(d) t_{pw} = 4.8ms');
    end


%% ------ F I G U R E  6 -------------------------------------------------
elseif flags.do_fig6

    if flags.do_plot
        [data,description] = data_wierstorf2013('fig6','plot');
    else
        [data,description] = data_wierstorf2013('fig6','noplot');
    end
    output.data = data;
    output.description = description;


%% ------ F I G U R E  7 -------------------------------------------------
elseif flags.do_fig7

    [data,description] = data_wierstorf2013('fig7','noplot');
    output.data = data;
    output.description = description;
    if flags.do_plot
        L = 2.85;
        conf.array = 'linear';
        conf.X0 = [0 0];
        conf.x0 = [];
        figure;
        subplot(1,3,1);
        conf.dx0 = L/2;
        quiver(data(:,1),data(:,2),data(:,3),data(:,4),25,'.b');
        hold on;
        draw_loudspeakers(secondary_source_positions(L,conf));
        hold on;
        plot(0,1,'*r');
        hold off;
        title(description{3,1});
        axis([-2.13 1.63 -2.2 1.2]);
        xlabel('x/m');
        ylabel('y/m');
        subplot(1,3,2);
        conf.dx0 = L/7;
        quiver(data(:,1),data(:,2),data(:,5),data(:,6),25,'.b');
        hold on;
        draw_loudspeakers(secondary_source_positions(L,conf));
        hold on;
        plot(0,1,'*r');
        hold off;
        title(description{5,1});
        axis([-2.13 1.63 -2.2 1.2]);
        xlabel('x/m');
        ylabel('y/m');
        subplot(1,3,3);
        conf.dx0 = L/14;
        quiver(data(:,1),data(:,2),data(:,7),data(:,8),25,'.b');
        hold on;
        draw_loudspeakers(secondary_source_positions(L,conf));
        hold on;
        plot(0,1,'*r');
        hold off;
        title(description{7,1});
        axis([-2.13 1.63 -2.2 1.2]);
        xlabel('x/m');
        ylabel('y/m');
    end


%% ------ F I G U R E  8 -------------------------------------------------
elseif flags.do_fig8

    s = [mfilename('fullpath'),'_fig8.mat'];

    if amtredofile(s,flags.redomode)
        % load HRTFs, see:
        % https://dev.qu.tu-berlin.de/projects/measurements/wiki/2010-11-kemar-anechoic
        [~,path] = download_hrtf('wierstorf2011_3m');
        load([path 'wierstorf2011_3m.mat']);
        % generate noise signal
        sig_noise = noise(44100/5,1,'white');
        % get only the -90 to 90 degree part of the hrtf set
        idx = (( irs.apparent_azimuth>-pi/2 & irs.apparent_azimuth<pi/2 & ...
            irs.apparent_elevation==0 ));
        irs = slice_irs(irs,idx);
        for ii=1:length(irs.apparent_azimuth)
            % generate noise coming from the given direction
            ir = get_ir(irs,[irs.apparent_azimuth(ii) 0 irs.distance]);
            sig = auralize_ir(ir,sig_noise);
            % calculate binaural parameters
            [fine, modulation, cfreqs, ild_tmp] = dietz2011(sig,44100);
            % unwrap ITD
            itd_tmp = ...
                dietz2011unwrapitd(fine.itd(:,1:12),ild_tmp(:,1:12),fine.f_inst,2.5);
            % calculate the mean about time of the binaural parameters and store
            % them
            itd(ii,:) = median(itd_tmp,1);
            phi(ii) = degree(irs.apparent_azimuth(ii));
        end
        save(save_format,s,'phi','itd');
    else
        load(s);
    end;

    output.phi = phi;
    output.itd = itd;

    if flags.do_plot
        % ------ Plotting ------
        figure;
        plot(phi,itd.*1000);
        axis([-90 0 0 0.9])
        xlabel('phi_{sound event}/deg');
        ylabel('interaural time difference/ms');
    end;

end;


%% ------ F I G U R E  9 -------------------------------------------------
if flags.do_fig9

    s = [mfilename('fullpath'),'_fig9.mat'];
      
    if amtredofile(s,flags.redomode)
        % load lookup table
        path = which('amtstart');
        lookup = ...
            load([path(1:end-10) 'modelstages/wierstorf2013itd2anglelookup.mat']);
        % load HRTFs, see:
        % https://dev.qu.tu-berlin.de/projects/measurements/wiki/2010-11-kemar-anechoic
        [~,path] = download_hrtf('wierstorf2011_3m');
        load([path 'wierstorf2011_3m.mat']);
        % generate noise signal
        sig_noise = noise(44100/5,1,'white');
        % get only the -90 to 90 degree part of the hrtf set
        idx = (( irs.apparent_azimuth>-pi/2 & irs.apparent_azimuth<pi/2 & ...
            irs.apparent_elevation==0 ));
        irs = slice_irs(irs,idx);
        for ii=1:length(irs.apparent_azimuth)
            % generate noise coming from the given direction
            ir = get_ir(irs,[irs.apparent_azimuth(ii) 0 irs.distance]);
            sig = auralize_ir(ir,sig_noise);
            phi_auditory_event(ii) = estimate_azimuth(sig,lookup,'dietz2011');
            phi_sound_event(ii) = degree(irs.apparent_azimuth(ii));
        end
        save(save_format,s,'phi_auditory_event','phi_sound_event');
    else
        load(s);
    end;

    output.phi_auditory_event = phi_auditory_event;
    output.phi_sound_event = phi_sound_event;

    if flags.do_plot
        % ------ Plotting ------
        figure;
        plot(phi_sound_event,phi_auditory_event-phi_sound_event);
        axis([-90 90 -3 3])
        xlabel('phi_{sound event}/deg');
        ylabel('phi_{auditory event}-phi_{sound event}/deg');
    end;
end


%% ------ F I G U R E  10 ------------------------------------------------
if flags.do_fig10

    s = [mfilename('fullpath'),'_fig10.mat'];

    % orientation of the listener (always to the front)
    phi = pi/2;
    % position of the virtual point source
    xs = [0 1];
    src = 'ps';
    % array size
    L = 2.85;

    % Y=1.5m
    X = -1.75:0.25:0;
    Y1 = 1.5;
    Y2 = 2.0;

    if amtredofile(s,flags.redomode)
        for ii=1:length(X)
            progressbar(ii,length(X));
            model_3_Y1(ii) = wierstorf2013(X,Y1,phi,xs,src,L,'wfs', ...
                                           'resolution',1, ...
                                           'nls',3, ...
                                           'array','linear', ...
                                           'showprogress',0);
            model_3_Y2(ii) = wierstorf2013(X,Y2,phi,xs,src,L,'wfs', ...
                                           'resolution',1, ...
                                           'nls',3, ...
                                           'array','linear', ...
                                           'showprogress',0);
            model_8_Y1(ii) = wierstorf2013(X,Y1,phi,xs,src,L,'wfs', ...
                                           'resolution',1, ...
                                           'nls',8, ...
                                           'array','linear', ...
                                           'showprogress',0);
            model_8_Y2(ii) = wierstorf2013(X,Y2,phi,xs,src,L,'wfs', ...
                                           'resolution',1, ...
                                           'nls',8, ...
                                           'array','linear', ...
                                           'showprogress',0);
            model_15_Y1(ii) = wierstorf2013(X,Y1,phi,xs,src,L,'wfs', ...
                                           'resolution',1, ...
                                           'nls',15, ...
                                           'array','linear', ...
                                           'showprogress',0);
            model_15_Y2(ii) = wierstorf2013(X,Y2,phi,xs,src,L,'wfs', ...
                                           'resolution',1, ...
                                           'nls',15, ...
                                           'array','linear', ...
                                           'showprogress',0);
        end
        save(save_format,s,'model_3_Y1','model_3_Y2','model_8_Y1','model_8_Y2','model_15_Y1','model_15_Y2');
    else
        load(s);
    end

    % get the human data
    [data,description] = data_wierstorf2013('fig10','noplot');
    
    output.model_3_Y1 = model_3_Y1;
    output.model_3_Y2 = model_3_Y2;
    output.model_8_Y1 = model_8_Y1;
    output.model_8_Y2 = model_8_Y2;
    output.model_15_Y1 = model_15_Y1;
    output.model_15_Y2 = model_15_Y2;
    output.data = data;
    output.description = description;

    if flags.do_plot

        data_wierstorf2013('fig10','plot');
        hold on;
        subplot(3,1,1)
        plot(data(:,1),model_3_Y1,'-b');
        plot(data(:,1),model_3_Y2,'-r');
        subplot(3,1,2)
        plot(data(:,1),model_8_Y1,'-b');
        plot(data(:,1),model_8_Y2,'-r');
        subplot(3,1,3)
        plot(data(:,1),model_15_Y1,'-b');
        plot(data(:,1),model_15_Y2,'-r');

    end
end


%% ------ F I G U R E  11a -----------------------------------------------
if flags.do_fig11a

    s = [mfilename('fullpath'),'_fig11a.mat'];
      
    % listening area
    X = [-2 2];
    Y = [-3.15 -0.15];
    % orientation of the listener (always to the front)
    phi = pi/2;
    % position of the virtual point source
    xs = [0 1];
    src = 'ps';
    % array size
    L = 2.85;
  
    if amtredofile(s,flags.redomode)
        fprintf(1,'\nWarning: this will take a long time!\n\n');
        % 3 speakers
        fprintf(1,'Calculating figure 1/6\n');
        [~,aud_event_3,~,xaxis_31,yaxis_31,x0_3] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                          'resolution',31, ...
                          'nls',3, ...
                          'array','linear');
        fprintf(1,'Calculating figure 2/6\n');
        [loc_error_3,~,~,xaxis_135,yaxis_135] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',3, ...
                'array','linear');
        % 8 speakers
        fprintf(1,'Calculating figure 3/6\n');
        [~,aud_event_8,~,~,~,x0_8] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',31, ...
                'nls',8, ...
                'array','linear');
        fprintf(1,'Calculating figure 4/6\n');
        loc_error_8 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',8, ...
                'array','linear');
        % 15 speakers
        fprintf(1,'Calculating figure 5/6\n');
        [~,aud_event_15,~,~,~,x0_15] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',31, ...
                'nls',15, ...
                'array','linear');
        fprintf(1,'Calculating figure 6/6\n');
        loc_error_15 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',15, ...
                'array','linear');
        save(save_format,s,'loc_error_3','loc_error_8','loc_error_15', ...
            'aud_event_3','aud_event_8','aud_event_15', ...
            'x0_3','x0_8','x0_15', ...
            'xaxis_31','yaxis_31','xaxis_135','yaxis_135');
    else
        load(s);
    end;

    output.loc_error_3 = loc_error_3;
    output.loc_error_8 = loc_error_8;
    output.loc_error_15 = loc_error_15;
    output.aud_event_3 = aud_event_3;
    output.aud_event_8 = aud_event_8;
    output.aud_event_15 = aud_event_15;
    output.x0_3 = x0_3;
    output.x0_8 = x0_8;
    output.x0_15 = x0_15;
    output.xaxis_31 = xaxis_31;
    output.yaxis_31 = yaxis_31;
    output.xaxis_135 = xaxis_135;
    output.yaxis_135 = yaxis_135;

    if flags.do_plot
        % ------ Plotting ------
        figure;
        subplot(2,3,1);
        [u,v,~] = pol2cart(rad(aud_event_3+90),ones(size(aud_event_3)), ...
            zeros(size(aud_event_3)));
        quiver(xaxis_31,yaxis_31,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_3,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,2);
        [u,v,~] = pol2cart(rad(aud_event_8+90),ones(size(aud_event_8)), ...
            zeros(size(aud_event_8)));
        quiver(xaxis_31,yaxis_31,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_8,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,3);
        [u,v,~] = pol2cart(rad(aud_event_15+90),ones(size(aud_event_15)), ...
            zeros(size(aud_event_15)));
        quiver(xaxis_31,yaxis_31,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_15,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,4)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_3'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_3,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,5)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_8'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_8,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,6)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_15'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_15,1);
        xlabel('x/m');
        ylabel('y/m');
    end;
end

%% ------ F I G U R E  11b -----------------------------------------------
if flags.do_fig11b

    s = [mfilename('fullpath'),'_fig11b.mat'];
      
    % listening area
    X = [-2 2];
    Y = [-3.15 -0.15];
    % orientation of the listener (always to the front)
    phi = pi/2;
    % position of the virtual point source
    xs = [0 -1];
    src = 'pw';
    % array size
    L = 2.85;
  
    if amtredofile(s,flags.redomode)
        fprintf(1,'\nWarning: this will take a long time!\n\n');
        % 3 speakers
        fprintf(1,'Calculating figure 1/6\n');
        [~,aud_event_3,~,xaxis_31,yaxis_31,x0_3] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                          'resolution',31, ...
                          'nls',3, ...
                          'array','linear');
        fprintf(1,'Calculating figure 2/6\n');
        [loc_error_3,~,~,xaxis_135,yaxis_135] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',3, ...
                'array','linear');
        % 8 speakers
        fprintf(1,'Calculating figure 3/6\n');
        [~,aud_event_8,~,~,~,x0_8] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',31, ...
                'nls',8, ...
                'array','linear');
        fprintf(1,'Calculating figure 4/6\n');
        loc_error_8 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',8, ...
                'array','linear');
        % 15 speakers
        fprintf(1,'Calculating figure 5/6\n');
        [~,aud_event_15,~,~,~,x0_15] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',31, ...
                'nls',15, ...
                'array','linear');
        fprintf(1,'Calculating figure 6/6\n');
        loc_error_15 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',15, ...
                'array','linear');
        save(save_format,s,'loc_error_3','loc_error_8','loc_error_15', ...
            'aud_event_3','aud_event_8','aud_event_15', ...
            'x0_3','x0_8','x0_15', ...
            'xaxis_31','yaxis_31','xaxis_135','yaxis_135');
    else
        load(s);
    end;

    output.loc_error_3 = loc_error_3;
    output.loc_error_8 = loc_error_8;
    output.loc_error_15 = loc_error_15;
    output.aud_event_3 = aud_event_3;
    output.aud_event_8 = aud_event_8;
    output.aud_event_15 = aud_event_15;
    output.x0_3 = x0_3;
    output.x0_8 = x0_8;
    output.x0_15 = x0_15;
    output.xaxis_31 = xaxis_31;
    output.yaxis_31 = yaxis_31;
    output.xaxis_135 = xaxis_135;
    output.yaxis_135 = yaxis_135;

    if flags.do_plot
        % ------ Plotting ------
        figure;
        subplot(2,3,1);
        [u,v,~] = pol2cart(rad(aud_event_3+90),ones(size(aud_event_3)), ...
            zeros(size(aud_event_3)));
        quiver(xaxis_31,yaxis_31,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_3,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,2);
        [u,v,~] = pol2cart(rad(aud_event_8+90),ones(size(aud_event_8)), ...
            zeros(size(aud_event_8)));
        quiver(xaxis_31,yaxis_31,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_8,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,3);
        [u,v,~] = pol2cart(rad(aud_event_15+90),ones(size(aud_event_15)), ...
            zeros(size(aud_event_15)));
        quiver(xaxis_31,yaxis_31,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_15,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,4)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_3'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_3,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,5)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_8'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_8,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,6)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_15'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_15,1);
        xlabel('x/m');
        ylabel('y/m');
    end;
end

%% ------ F I G U R E  12a -----------------------------------------------
if flags.do_fig12a

    s = [mfilename('fullpath'),'_fig12a.mat'];
      
    % listening area
    X = [-2.1 2.1];
    Y = [-2.1 2.1];
    % orientation of the listener (always to the front)
    phi = pi/2;
    % position of the virtual point source
    xs = [0 2.5];
    src = 'ps';
    % array size
    L = 3;
  
    if amtredofile(s,flags.redomode)
        fprintf(1,'\nWarning: this will take a long time!\n\n');
        % 14 speakers
        fprintf(1,'Calculating figure 1/6\n');
        [~,aud_event_14,~,xaxis_21,yaxis_21,x0_14] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',21, ...
                'nls',14, ...
                'array','circle');
        fprintf(1,'Calculating figure 2/6\n');
        [loc_error_14,~,~,xaxis_135,yaxis_135] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',14, ...
                'array','circle');
        % 28 speakers
        fprintf(1,'Calculating figure 3/6\n');
        [~,aud_event_28,~,~,~,x0_28] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',21, ...
                'nls',28, ...
                'array','circle');
        fprintf(1,'Calculating figure 4/6\n');
        loc_error_28 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',28, ...
                'array','circle');
        % 56 speakers
        fprintf(1,'Calculating figure 5/6\n');
        [~,aud_event_56,~,~,~,x0_56] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',21, ...
                'nls',56, ...
                'array','circle');
        fprintf(1,'Calculating figure 6/6\n');
        loc_error_56 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',56, ...
                'array','circle');
        save(save_format,s,'loc_error_14','loc_error_28','loc_error_56', ...
            'aud_event_14','aud_event_28','aud_event_56', ...
            'x0_14','x0_28','x0_56', ...
            'xaxis_21','yaxis_21','xaxis_135','yaxis_135');
    else
        load(s);
    end;

    output.loc_error_14 = loc_error_14;
    output.loc_error_28 = loc_error_28;
    output.loc_error_56 = loc_error_56;
    output.aud_event_14 = aud_event_14;
    output.aud_event_28 = aud_event_28;
    output.aud_event_56 = aud_event_56;
    output.x0_14 = x0_14;
    output.x0_28 = x0_28;
    output.x0_56 = x0_56;
    output.xaxis_21 = xaxis_21;
    output.yaxis_21 = yaxis_21;
    output.xaxis_135 = xaxis_135;
    output.yaxis_135 = yaxis_135;

    if flags.do_plot
        % ------ Plotting ------
        figure;
        subplot(2,3,1);
        [u,v,~] = pol2cart(rad(aud_event_14+90),ones(size(aud_event_14)), ...
            zeros(size(aud_event_14)));
        quiver(xaxis_21,yaxis_21,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_14,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,2);
        [u,v,~] = pol2cart(rad(aud_event_28+90),ones(size(aud_event_28)), ...
            zeros(size(aud_event_28)));
        quiver(xaxis_21,yaxis_21,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_28,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,3);
        [u,v,~] = pol2cart(rad(aud_event_56+90),ones(size(aud_event_56)), ...
            zeros(size(aud_event_56)));
        quiver(xaxis_21,yaxis_21,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_56,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,4)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_14'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_14,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,5)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_28'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_28,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,6)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_56'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_56,1);
        xlabel('x/m');
        ylabel('y/m');
    end;
end

%% ------ F I G U R E  12b -----------------------------------------------
if flags.do_fig12b

    s = [mfilename('fullpath'),'_fig12b.mat'];
      
    % listening area
    X = [-2.1 2.1];
    Y = [-2.1 2.1];
    % orientation of the listener (always to the front)
    phi = pi/2;
    % position of the virtual point source
    xs = [0 -1];
    src = 'pw';
    % array size
    L = 3;
  
    if amtredofile(s,flags.redomode)
        fprintf(1,'\nWarning: this will take a long time!\n\n');
        % 14 speakers
        fprintf(1,'Calculating figure 1/6\n');
        [~,aud_event_14,~,xaxis_21,yaxis_21,x0_14] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',21, ...
                'nls',14, ...
                'array','circle');
        fprintf(1,'Calculating figure 2/6\n');
        [loc_error_14,~,~,xaxis_135,yaxis_135] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',14, ...
                'array','circle');
        % 28 speakers
        fprintf(1,'Calculating figure 3/6\n');
        [~,aud_event_28,~,~,~,x0_28] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',21, ...
                'nls',28, ...
                'array','circle');
        fprintf(1,'Calculating figure 4/6\n');
        loc_error_28 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',28, ...
                'array','circle');
        % 56 speakers
        fprintf(1,'Calculating figure 5/6\n');
        [~,aud_event_56,~,~,~,x0_56] = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',21, ...
                'nls',56, ...
                'array','circle');
        fprintf(1,'Calculating figure 6/6\n');
        loc_error_56 = ...
            wierstorf2013(X,Y,phi,xs,src,L,'wfs', ...
                'resolution',135, ...
                'nls',56, ...
                'array','circle');
        save(save_format,s,'loc_error_14','loc_error_28','loc_error_56', ...
            'aud_event_14','aud_event_28','aud_event_56', ...
            'x0_14','x0_28','x0_56', ...
            'xaxis_21','yaxis_21','xaxis_135','yaxis_135');
    else
        load(s);
    end;

    output.loc_error_14 = loc_error_14;
    output.loc_error_28 = loc_error_28;
    output.loc_error_56 = loc_error_56;
    output.aud_event_14 = aud_event_14;
    output.aud_event_28 = aud_event_28;
    output.aud_event_56 = aud_event_56;
    output.x0_14 = x0_14;
    output.x0_28 = x0_28;
    output.x0_56 = x0_56;
    output.xaxis_21 = xaxis_21;
    output.yaxis_21 = yaxis_21;
    output.xaxis_135 = xaxis_135;
    output.yaxis_135 = yaxis_135;

    if flags.do_plot
        % ------ Plotting ------
        figure;
        subplot(2,3,1);
        [u,v,~] = pol2cart(rad(aud_event_14+90),ones(size(aud_event_14)), ...
            zeros(size(aud_event_14)));
        quiver(xaxis_21,yaxis_21,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_14,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,2);
        [u,v,~] = pol2cart(rad(aud_event_28+90),ones(size(aud_event_28)), ...
            zeros(size(aud_event_28)));
        quiver(xaxis_21,yaxis_21,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_28,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,3);
        [u,v,~] = pol2cart(rad(aud_event_56+90),ones(size(aud_event_56)), ...
            zeros(size(aud_event_56)));
        quiver(xaxis_21,yaxis_21,u',v',0.5);
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_56,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,4)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_14'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_14,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,5)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_28'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_28,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,6)
        imagesc(xaxis_135,yaxis_135,abs(loc_error_56'));
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_56,1);
        xlabel('x/m');
        ylabel('y/m');
    end;
end
