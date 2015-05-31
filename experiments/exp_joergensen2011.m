function exp_joergensen2011(varargin)
%EXP_JOERGENSEN2011 Figures from Jørgensen and Dau (2011)
%   Usage: output = exp_joergensen2011(flag)
%
%   `exp_joergensen2011(flag)` reproduces the results for the figure given
%   by flag from the Jørgensen and Dau (2011) paper. 
%    
%   The following flags can be specified;
%
%     'plot'     Plot the specified figure from Jørgensen and Dau (2011). This is
%                the default. 
%
%     'noplot'   Don't plot, only return data.
%
%
%     'redo'     Always recalculate the experiment results.
%
%     'cached'   Always use the cached version. Default.
%
%     'fig5'     Plot Fig. 5 (Jørgensen and Dau, 2011).
%
%     'fig6'     Plot Fig. 6 (Jørgensen and Dau, 2011). 
% 
%
%   Examples:
%   ---------
%
%   To display Figure 5 use :::
%
%     exp_joergensen2011('fig5');
%
%   To display Figure 6 use :::
%
%     exp_joergensen2011('fig6');
%
%   ---------
%
%   Please cite Jørgensen and Dau (2011) if you use
%   this model.
%
%   See also: joergensen2011, plotjoergensen2011, exp_joergensen2011
%
%   References: joergensen2011predicting

definput.import={'amtcache'};
definput.flags.type={'missingflag','fig5','fig6'};
definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

NSpeechsamples = 10; % specify the number of speech samples to be used fro the simulations. The simulation takes longer the more samples are used. A minimum of 50 should be used for final validation.  

%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5;

  dSRT = amtcache('get', ['fig5_' num2str(NSpeechsamples) 'sntcs'], flags.cachemode);

  if isempty(dSRT)
    dSRT = joergensen2011sim(NSpeechsamples,'fig5');
    amtcache('set',['fig5_' num2str(NSpeechsamples) 'sntcs'],dSRT);
  end;
        
  if flags.do_plot
    plotjoergensen2011(dSRT,'fig5');
  end
end;

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6;

  dSRT = amtcache('get', ['fig6_' num2str(NSpeechsamples) 'sntcs'], flags.cachemode);
  if isempty(dSRT)

    dSRT = joergensen2011sim(NSpeechsamples,'fig6');
		amtcache('set',['fig6_' num2str(NSpeechsamples) 'sntcs'],dSRT);
  end;
        
  if flags.do_plot
    plotjoergensen2011(dSRT,'fig6');
  end
end;
