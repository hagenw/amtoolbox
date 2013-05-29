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
if !which('SFS_start')
    error(['%s: you need to install the Sound-Field-Synthesis Toolbox.\n', ...
        'You can download it at https://github.com/sfstoolbox/sfs.\n', ...
        'The results in the paper are verified up to revision a8914700a4'], ...
        upper(mfilename));
end
conf = SFS_config;

%% ------ FIG 1 -----------------------------------------------------------
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
        s = load(s);
        loc_error = s.loc_error;
        aud_event = s.aud_event;
        xaxis = s.xaxis;
        yaxis = s.yaxis;
        x0 = s.x0;
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

elseif flags.do_fig3

    s = [mfilename('fullpath'),'_fig3.mat'];

    conf = SFS_config;  
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
end;
