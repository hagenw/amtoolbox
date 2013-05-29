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
%     'fig1'  Reproduce Fig.1 from Wierstorf (2013). The localization error
%             for a typical stereophony setup is calculated and shown for the
%             whole listening are, sampled with 21x21 point.
%
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
definput.flags.type={'missingflag','fig1','fig3','fig6','fig7','fig8',...
                    'fig9','fig10','fig11','fig12'};

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
        'The results in the paper are verified up to revision a8914700a4'], ...
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
            wierstorf2013(X,Y,phi,xs,src,L,res,'wfs', ...
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
        imagesc(xaxis_135,yaxis_135,loc_error_3');
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_3,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,5)
        imagesc(xaxis_135,yaxis_135,loc_error_8');
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_8,1);
        xlabel('x/m');
        ylabel('y/m');
        subplot(2,3,6)
        imagesc(xaxis_135,yaxis_135,loc_error_15');
        turn_imagesc;
        axis([-2.13 2.13 -3.3 0.2])
        draw_loudspeakers(x0_15,1);
        xlabel('x/m');
        ylabel('y/m');
    end;
end
