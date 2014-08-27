function exp_joergensen2011(NSpeechsamples,varargin)
%EXP_JOERGENSEN2011 Figures from Jørgensen and Dau (2011)
%   Usage: output = exp_joergensen2011(flag)
%
%   `exp_joergensen2011(flag)` reproduces the results for the figure given
%   by flag from the Jørgensen and Dau (2011) paper. 
% 
%   Inputs:
%      NSpeechsamples: specify the number of speech samples to be used fro the simulations. The simulation takes longer the more samples are used. A minimum of 50 should be used for final validation.  
%    
%   The following flags can be specified;
%
%     'plot'     Plot the specified figure from Jørgensen and Dau (2011). This is
%                the default. 
%
%     'noplot'   Don't plot, only return data.
%
%
%     'auto '    Redo the experiment if a cached dataset does not exist. This
%                is the default.
% 
%     'refresh'  Always recalculate the experiment.
%
%     'cached'   Always use the cached version. This throws an error if the
%                file does not exist.
%
%     'fig5'     Plot Fig. 5 (Jørgensen and Dau, 2011).
%
%     'fig6'     Plot Fig. 6 (Jørgensen and Dau, 2011). 
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

definput.import={'amtredofile'};
definput.flags.type={'fig5','fig6'};
definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

save_format='-v6';
%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5;

          

  s = [mfilename('fullpath'),'_fig5_' num2str(NSpeechsamples)  'sntcs.mat'];

  if amtredofile(s,flags.redomode)

    [dSRT] = joergensen2011sim(NSpeechsamples,'fig5');

    save(s,'dSRT',save_format);
  else

    s = load(s);
   dSRT = s.dSRT;
  end;
        
  if flags.do_plot
    plotjoergensen2011(dSRT,'fig5');
  end
end;

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6;

          

  s = [mfilename('fullpath'),'_fig6_' num2str(NSpeechsamples)  'sntcs.mat'];

  if amtredofile(s,flags.redomode)

    [dSRT] = joergensen2011sim(NSpeechsamples,'fig6');

    save(s,'dSRT',save_format);
  else

    s = load(s);
   dSRT = s.ans;
  end;
        
  if flags.do_plot
    plotjoergensen2011(dSRT,'fig6');
  end
end;
