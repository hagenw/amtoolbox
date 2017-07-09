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
%     'redo'     Recalculate the experiment results.
%
%     'cached'   Use cached results. Default.
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
%
%   ---------
%
%   Please cite Rønne et al. (2012) and Zilany and Bruce (2007) if you use
%   this model.
%
definput.import={'amt_cache'};
definput.flags.type={'missingflag','fig5','fig6','fig7'};
definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
		flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}), ...
				sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
		error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5
  
  waveVamp = 0;

  stim_level = 40:10:100; % Default stimulus levels

  [waveVlat, click_latency] = amt_cache('get','fig5',flags.cachemode);  
  if isempty(waveVlat);
    [click_amplitude, click_latency]    = roenne2012click(stim_level);     
    waveVlat = roenne2012tonebursts(stim_level);    
    amt_cache('set','fig5',waveVlat,click_latency);
  end;
  
  if flags.do_plot;
    plot_roenne2012tonebursts(waveVlat,click_latency);
  end  ;  

end;

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6;

  stim_level    = (20:20:60)+35.2;
  
  % Default chirp numbers. 1 = click, 2 to 6 = chirp 1 to 5.
  chirp_number  = 1:6;              

  [waveVamp, waveVlat] = amt_cache('get','fig6',flags.cachemode);
  if isempty(waveVamp)
    [waveVamp, waveVlat] = roenne2012chirp(stim_level, chirp_number);
    amt_cache('set','fig6',waveVamp,waveVlat);
  end;
        
  if flags.do_plot
    plot_roenne2012chirp(waveVamp, waveVlat,'amponly');
  end
end;

%% ------ FIG 7 -----------------------------------------------------------
if flags.do_fig7;

  stim_level    = (20:20:60)+35.2;
  
  % Default chirp numbers. 1 = click, 2 to 6 = chirp 1 to 5.
  chirp_number  = 1:6;              

  [waveVamp, waveVlat] = amt_cache('get','fig7',flags.cachemode);
  if isempty(waveVamp)
    [waveVamp, waveVlat] = roenne2012chirp(stim_level, chirp_number);
    amt_cache('set','fig7',waveVamp,waveVlat);
  end;
        
  if flags.do_plot
    plot_roenne2012chirp(waveVamp, waveVlat,'latonly');
  end
end;

