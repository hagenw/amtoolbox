function exp_joergensen2013(NSpeechsamples, varargin)
%EXP_JOERGENSEN2013 Figures from Jørgensen, Ewert and Dau (2013)
%   Usage: output = exp_joergensen2013(flag)
%
%   `exp_joergensen2013(NSpeechsamples, flag)` reproduces the results for
%   figure 2 given by `flag` from Jørgensen, Ewert and Dau (2013). The input
%   *NSpeechsamples* specifies the number of speech samples to be used for the
%   simulations. The simulation takes longer the more samples are used. A
%   minimum of 50 should be used for final validation.
% 
%   The following flags can be specified;
%
%     'plot'     Plot the specified figure from Jørgensen, Ewert and Dau (2013). This is
%                the default.
%
%     'noplot'   Don't plot, only return data.
%
%     'redo'     Always recalculate the experiment.
%
%     'cached'   Always use the cached version. Default.
% 
%   Only one of the following:
%     'fig2_simAll'     Simualte all data in fig2 of Jørgensen, Ewert and Dau (2013).
%
%     'fig2_specsub'    Simulate and plot the conditions with spectral subtraction shown in fig2 (Jørgensen, Ewert and Dau, 2013) 
% 
%     'fig2_reverb'     Simulate and plot the conditions with reverberation shown in fig2 (Jørgensen, Ewert and Dau, 2013) 
% 
%     'fig2_kjems2009'  Simulate and plot the data from Kjems et al (2009) shown in fig2 (Jørgensen, Ewert and Dau, 2013) 
% 
%     'fig2_FP1990'     Simulate and plot the data from Festen and Plomp (1990) shown in fig2 (Jørgensen, Ewert and Dau, 2013) 
% 
%     'fig2_Jetal2013'  Simulate and plot the new data shown in fig2 (Jørgensen, Ewert and Dau, 2013) 
% 
%     'fig2_OrgSim'     plot the simulation shown in fig2 (Jørgensen, Ewert and Dau, 2013) 
%
%
%   Examples:
%   ---------
%
%   To simulate all conditions and display Figure 2 use :::
%
%     exp_joergensen2013('fig2_simAll');
%  
%   See also: joergensen2013, joergensen2013sim, plotjoergensen2013
%
%   ---------
%
%   Please cite Jørgensen et al. (2013) if you use
%   this model.
% 
%   References: joergensen2013 joergensen2011predicting kjems2009role festen1990effects

definput.import={'amtcache'};
definput.flags.type={'missingflag','fig2_simAll','fig2_specsub','fig2_reverb','fig2_kjems2009','fig2_FP1990','fig2_OrgSim','fig2_Jetal2013'};
definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

%% ------ FIG 2 -----------------------------------------------------------
if flags.do_fig2_simAll;
        
    SRTs = amtcache('get',['fig2_' num2str(NSpeechsamples)  'sntcs'], flags.cachemode);
    if isempty(SRTs)
                       
        [SRTs_specsub specsubConds] = joergensen2013sim(NSpeechsamples,'JandD2011specsub');
        [SRTs_reverb reverbConds] = joergensen2013sim(NSpeechsamples,'JandD2011reverb');
        [SRTs_simJetal2013 Jetal2013Conds] = joergensen2013sim(NSpeechsamples,'Jetal2013');
        [SRTs_simFP1990 FP1990Conds] = joergensen2013sim(NSpeechsamples,'FP1990');
        [SRTs_simKjems2009 Kjems2009Conds] = joergensen2013sim(NSpeechsamples,'Kjems2009');
        
        SRTs.simSRTs_Jetal2013 = SRTs_simJetal2013;
        SRTs.simSRTs_reverb = SRTs_reverb;
        SRTs.simSRTs_specsub = SRTs_specsub;
        SRTs.simSRTs_Kjems2009 = SRTs_simKjems2009;
        SRTs.simSRTs_FP1990 = SRTs_simFP1990;
        
        amtcache('set',['fig2_' num2str(NSpeechsamples)  'sntcs'],SRTs);
    end;
    
    if flags.do_plot
        plotjoergensen2013(SRTs,'fig2');
    end
end;

if flags.do_fig2_specsub;
    
    SRTs = amtcache('get',['fig2_specsub_' num2str(NSpeechsamples)  'sntcs'], flags.cachemode);
    if isempty(SRTs)
               
        [SRTs_specsub specsubConds] = joergensen2013sim(NSpeechsamples,'JandD2011specsub');
        SRTs.simSRTs_specsub = SRTs_specsub;
        SRTs.simSRTs_Kjems2009 = [ NaN NaN NaN NaN ] ;
        SRTs.simSRTs_FP1990 =[ NaN NaN NaN ];
        SRTs.simSRTs_Jetal2013 = [ NaN NaN NaN  ];
        SRTs.simSRTs_reverb =[ NaN NaN NaN NaN NaN  ];

        amtcache('set', ['fig2_specsub_' num2str(NSpeechsamples)  'sntcs'], SRTs);
    end;
    
    if flags.do_plot
        plotjoergensen2013(SRTs,'fig2');
    end
end;

if flags.do_fig2_reverb;
 
    SRTs = amtcache('get',['fig2_reverb_' num2str(NSpeechsamples)  'sntcs'], flags.cachemode);
    if isempty(SRTs)
                       
        [SRTs_reverb reverbConds] = joergensen2013sim(NSpeechsamples,'JandD2011reverb');
        SRTs.simSRTs_specsub = [ NaN NaN NaN NaN NaN NaN ];
        SRTs.simSRTs_Kjems2009 = [ NaN NaN NaN NaN ];
        SRTs.simSRTs_FP1990 = [ NaN NaN NaN ];
        SRTs.simSRTs_Jetal2013 = [ NaN NaN NaN ];
        SRTs.simSRTs_reverb = SRTs_reverb;

        amtcache('set', ['fig2_reverb_' num2str(NSpeechsamples)  'sntcs'], SRTs);
    end;
    
    if flags.do_plot
        plotjoergensen2013(SRTs,'fig2');
    end
end;

if flags.do_fig2_kjems2009;
   
    SRTs = amtcache('get',['fig2_Kjems2009_' num2str(NSpeechsamples)  'sntcs'], flags.cachemode);
    if isempty(SRTs)
                       
        [SRTs_Kjems2009 Kjems2009Conds] = joergensen2013sim(NSpeechsamples,'Kjems2009');
         SRTs.simSRTs_specsub = [ NaN NaN NaN NaN NaN NaN ];
         SRTs.simSRTs_Kjems2009 = SRTs_Kjems2009;
        SRTs.simSRTs_FP1990 = [ NaN NaN NaN  ];
        SRTs.simSRTs_Jetal2013 = [ NaN NaN NaN  ];
        SRTs.simSRTs_reverb = [ NaN NaN NaN NaN NaN  ];

        amtcache('set', ['fig2_Kjems2009_' num2str(NSpeechsamples)  'sntcs'], SRTs);
    end;
    
    if flags.do_plot
        plotjoergensen2013(SRTs,'fig2');
    end
end;

if flags.do_fig2_FP1990;
   
    SRTs = amtcache('get',['fig2_FP1990_' num2str(NSpeechsamples)  'sntcs'], flags.cachemode);
    if isempty(SRTs)
                      
        [SRTs_FP1990 FP1990Conds] = joergensen2013sim(NSpeechsamples,'FP1990');
        SRTs.simSRTs_specsub = [ NaN NaN NaN NaN NaN NaN ];
        SRTs.simSRTs_Kjems2009 = [ NaN NaN NaN NaN  ];
        SRTs.simSRTs_FP1990 = SRTs_FP1990;
        SRTs.simSRTs_Jetal2013 = [ NaN NaN NaN  ];
        SRTs.simSRTs_reverb = [ NaN NaN NaN NaN NaN  ];

        amtcache('set', ['fig2_FP1990_' num2str(NSpeechsamples)  'sntcs'], SRTs);
    end;
    
    if flags.do_plot
        plotjoergensen2013(SRTs,'fig2');
    end
end;

if flags.do_fig2_Jetal2013;
   
    SRTs = amtcache('get',['fig2_Jetal2013_' num2str(NSpeechsamples)  'sntcs'], flags.cachemode);
    if isempty(SRTs)
                      
        [SRTs_Jetal2013 Jetal2013Conds] = joergensen2013sim(NSpeechsamples,'Jetal2013');
        SRTs.simSRTs_specsub = [ NaN NaN NaN NaN NaN NaN ];
        SRTs.simSRTs_Kjems2009 = [ NaN NaN NaN NaN  ];
        SRTs.simSRTs_FP1990 = [ NaN NaN NaN  ];
        SRTs.simSRTs_Jetal2013 = SRTs_Jetal2013;
        SRTs.simSRTs_reverb = [ NaN NaN NaN NaN NaN  ];

        amtcache('set', ['fig2_Jetal2013_' num2str(NSpeechsamples)  'sntcs'], SRTs);
    end;
    
    if flags.do_plot
        plotjoergensen2013(SRTs,'fig2');
    end
end;

if  flags.do_fig2_OrgSim
    %load('plotting_jasa2012_final_predictionsJetal2013_fig2')
    %if flags.do_plot
    %    plotjoergensen2013(SRTs,'fig2');
    %end
end


