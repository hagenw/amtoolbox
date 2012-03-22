function [waveVamp, waveVlat] = exp_roenne2012(varargin)
%EXP_ROENNE2012 Figures from Rønne et al. (2012)
%   Usage: output = exp_roenne2012(flag)
%
%   `exp_roenne2012(flag)` reproduces the results for the figure given
%   by flag from the Rønne et al. (2012) paper. Outputs are the ABR wave V
%   amplitude and latency of all datapoints in that given figure.
%   
%   The following flags can be specified;
%
%     'plot'     Plot the specified figure from Rønne et al. (2012). This is
%                the default. 
%
%     'noplot'   Don't plot, only return data.
%
%     'plot2'    Plot extra figures for all individual simulated points.
%                Note that this creates lots of extra figures (3 for each
%                simulated data point)
%
%     'auto '    Re-calculate the file if it does not exist. Return 1 if the
%                file exist, otherwise 0. This is the default
%
%     'refresh'  Always recalculate the file.
%
%     'cached'   Always use the cached version. Throws an error if the
%                file does not exist.
%
%     'fig5'     Plot Fig. 5 (Rønne et al., 2012). Latency of simulated ABR
%                wave V's compared to Neely et al. (1988) and Harte et al.
%                (2009) reference data.
%
%     'fig6'     Plot Fig. 6 (Rønne et al., 2012). Amplitude of simulated
%                wave V compared to Elberling et al. (2010) reference data.
%
%     'fig7'     Plot Fig. 7 (Rønne et al., 2012). Latency of simulated wave
%                V compared to Elberling et al. (2010) reference data.
%
%   Examples:
%   ---------
%
%   To display Figure 5 use :::
%
%     exp_roenne2012('fig5');
%
%   To display Figure 6 use :::
%
%     exp_roenne2012('fig6');
%
%   To display Figure 7 use :::
%
%     exp_roenne2012('fig7');
%
%   References: roenne2012modeling elberling2010evaluating neely1988latency harte2009comparison

definput.import={'amtredofile'};
definput.flags.type={'fig5','fig6','fig7'};
definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

save_format='-v6';

%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5
  
  s = [mfilename('fullpath'),'_fig5.mat'];

  waveVamp = 0;

  stim_level                          = 40:10:100; % Default stimulus levels
  
  if amtredofile(s,flags.redomode)

    [click_amplitude, click_latency]    = roenne2012_click(stim_level);     

    waveVlat = roenne2012_tonebursts(stim_level);
    
    save(s,'waveVlat','click_latency',save_format);
  
  else
    
    s = load(s);
    waveVlat      = s.waveVlat;
    click_latency = s.click_latency;

  end;
  
  if flags.do_plot
    plotroenne2012_tonebursts(waveVlat,click_latency)
  end    

end;

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6

  stim_level    = (20:20:60)+35.2;
  
  % Default chirp numbers. 1 = click, 2 to 6 = chirp 1 to 5.
  chirp_number  = 1:6;              

  s = [mfilename('fullpath'),'_fig6.mat'];

  if amtredofile(s,flags.redomode)

    [waveVamp, waveVlat] = roenne2012_chirp(stim_level, chirp_number);

    save(s,'waveVamp','waveVlat',save_format);
  else

    s = load(s);
    waveVamp = s.waveVamp;
    waveVlat = s.waveVlat;
  end;
        
  if flags.do_plot
    plotroenne2012_chirp(waveVamp, waveVlat, flags.type)
  end
end;

if flags.do_fig7

  stim_level    = (20:20:60)+35.2;
  
  % Default chirp numbers. 1 = click, 2 to 6 = chirp 1 to 5.
  chirp_number  = 1:6;              

  s = [mfilename('fullpath'),'_fig7.mat'];

  if amtredofile(s,flags.redomode)

    [waveVamp, waveVlat] = roenne2012_chirp(stim_level, chirp_number);

    save(s,'waveVamp','waveVlat',save_format);
  else

    s = load(s);
    waveVamp = s.waveVamp;
    waveVlat = s.waveVlat;
  end;
        
  if flags.do_plot
    plotroenne2012_chirp(waveVamp, waveVlat, flags.type)
  end
end;

