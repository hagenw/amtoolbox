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
      
    % listener positions
    X = [-2 2];
    Y = [0.15 -3.15];
    conf.xyresolution = 21;
    [xx,yy,x,y] = xy_grid(X,Y,conf);
    % orientation of the listener (always to the front)
    phi = -pi/2;
    % position of the virtual point source
    xs = [0 0];
    src = 'ps';
    % distance between the loudspeaker ion the stereo setup
    L = 2;
    method = 'stereo';
  
    if amtredofile(s,flags.redomode)

        % estimate the localization
        for ii=1:length(x)
            [loc_error(ii),aud_event(ii),sound_event(ii),x0] = ...
                wierstorf2013([xx(ii) yy(ii)],phi,xs,src,L,'stereo');
        end

        save(s,'loc_error','aud_event','x0',save_format);
    else
        s = load(s);
        loc_error = s.loc_error;
        aud_event = s.aud_event;
        x0 = s.x0;
    end;
    
    if flags.do_plot
        % ------ Plotting ------
        figure;
        imagesc(x,y,degree(loc_error));
        draw_loudspeakers(x0,[1 1]);
        xlabel('x/m');
        ylabel('y/m');
        title('sweet spot in stereophony');
    end;

end;
