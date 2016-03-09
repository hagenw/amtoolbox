function varargout = exp_baumgartner2016(varargin)
%EXP_BAUMGARTNER2016 evaluation of baumgartner2016 model
%   Usage: data = exp_baumgartner2016(flag)
%
%   `exp_baumgartner2016(flag)` reproduces figures of the study from 
%   Baumgartner et al. (2014).
%
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'autorefresh'  Re-calculate the file if it does not exist. Return 1 if the
%                file exist, otherwise 0. This is the default
%
%     'refresh'  Always recalculate the file.
%
%     'cached'   Always use the cached version. Throws an error if the
%                file does not exist.
%
%     'baseline_ex' 
%               Prediction examples. Actual responses and response predictions 
%               for three exemplary listeners when listening to median-plane 
%               targets in the baseline condition. 
%               Actual response angles are shown as open circles. 
%               Probabilistic response predictions are encoded by brightness 
%               according to the color bar to the right. Actual (A:) and 
%               predicted (P:) quadrant error rates (QE) and local polar 
%               RMS errors (PE) are listed above each panel.
%               
%     'parameterization'
%               Model parametrization. 
%
%     'spatstrat_ex' 
%               Effect of band limitation and spectral warping. Actual 
%               responses and response predictions for listener NH12 when 
%               listening to broadband (BB), low-pass filtered (LP), or 
%               spectrally warped (W) DTFs of the median plane. Data were 
%               pooled within $\pm15^\circ$ of lateral angle.
%
%     'spatstrat' 
%               Effect of band limitation and spectral warping. 
%               Listeners were tested with broadband (BB), low-pass 
%               filtered (LP), and spectrally warped (W) DTFs. 
%               Actual: experimental results from majdak2013. 
%               Part.: Model predictions for the actual eight participants 
%               based on the actually tested target positions. Pool: Model 
%               predictions for our listener pool based on all possible 
%               target positions. Symbols and whiskers show median values 
%               and inter-quartile ranges, respectively. Symbols were 
%               horizontally shifted to avoid overlaps. Dotted horizontal 
%               lines represent chance rate. Correlation coefficients, $r$, 
%               and prediction residues, $e$, specify the correspondence 
%               between actual and predicted listener-specific performances.
%
%     'numchan_ex'
%               Effect of spectral resolution in terms of varying the number 
%               of spectral channels of a channel vocoder. Actual responses 
%               and response predictions for exemplary listener NH12. 
%               Results for 24, 9, and 3 channels are shown. All other 
%               conventions are as in Fig.3.
%
%     'numchan' 
%               Effect of spectral resolution in terms of varying the number 
%               of spectral channels of a channel vocoder. Actual experimental 
%               results are from goupell2010numchan. Stimulation with broadband 
%               click trains (CL) represents an unlimited number of channels. 
%               All other conventions are as in Fig.6. 
%
%     'baseline' 
%               Baseline performance as a function of the magnitude of the 
%               lateral response angle. Symbols and whiskers show median 
%               values and inter-quartile ranges, respectively. Open symbols 
%               represent actual and closed symbols predicted results. Symbols 
%               were horizontally shifted to avoid overlaps. Triangles with 
%               dashed lines show predictions (P) of the model without the 
%               sensomotoric mapping (SMM) stage.  
%
%     'hearingthreshold' 
%               Estimation of hearing thresholds in correspondence to OHC dysfunctions.  
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2014
%
%   3) Statistics Toolbox for Matlab (for some of the figures)
%
%   Examples:
%   ---------
%
%   To display the baseline prediction examples use :::
%
%     exp_baumgartner2016('baseline_ex');
%
%   To display the baseline prediction use :::
%
%     exp_baumgartner2016('baseline');
%
%   To display parametrization results use :::
%
%     exp_baumgartner2016('parametrization');
%
%   To display spatstrat prediction examples use :::
%
%     exp_baumgartner2016('spatstrat_ex');
%
%   To display spatstrat prediction use :::
%
%     exp_baumgartner2016('spatstrat');
%
%   To display numchan prediction examples use :::
%
%     exp_baumgartner2016('numchan_ex');
%
%   To display numchan prediction use :::
%
%     exp_baumgartner2016('numchan');
%
%   See also: baumgartner2014 data_baumgartner2016
%
%   References: baumgartner2014modeling majdak2013spatstrat goupell2010numchan  
%   


% AUTHOR: Robert Baumgartner


%% ------ Check input options --------------------------------------------

definput.import={'amtcache','localizationerror','baumgartner2014pmv2ppp'};

definput.flags.experiment = {'missingflag','baseline','baseline_ex',...
   'baseline_lat','parametrization',...
   'spatstrat_ex','spatstrat','numchan_ex','numchan',...
   'sabin2005','impairment','effectOnCues','dynrangecheck',...
   'localevel','hearingthreshold'};
     
definput.keyvals.ModelSettings = {};
 
% Figure Settings
definput.flags.plot = {'plot','noplot'};
definput.keyvals.FontSize = 10;
definput.keyvals.MarkerSize = 6;
definput.keyvals.gap = 0;
definput.keyvals.marg_h = [.08,.05];
definput.keyvals.marg_w = .05;
definput.keyvals.TickLength = [0.02,0.04];

% For: sabin2005
definput.keyvals.SL2SPL = 10;
definput.flags.gainInDeg = {'','gainInDeg'};

% For: impairment
definput.keyvals.subjects = [];
definput.flags.impairment={'FTlabel','noFTlabel'};

% For: localevel
definput.flags.TL = {'','TLisSPL'};
definput.flags.localevelplot = {'performance','pmv'};

% For: parametrization
definput.flags.parametrize = {'gamma','mrsandgamma'};

% For: dynrangecheck
definput.flags.dynrangecheck = {'','dynrangeDiff','dprime'};
definput.flags.dynrangecheck_comb = {'separate','combined'};

definput.flags.effectOnCues={'OHC','FT'};

% General: availability of Statistics Toolbox
definput.flags.statistics = {'stat','nostat'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

model.definput.import={'baumgartner2016'};
[model.flags,model.kv] = ltfatarghelper({},model.definput,kv.ModelSettings);

errorflag = [flags.errorflag,flags.ppp];

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.experiment{2:end-2}),...
             sprintf('%s or %s',definput.flags.experiment{end-1},definput.flags.experiment{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% Define cache name according to settings for auditory periphery model

cachename = ['g' num2str(model.kv.gamma,'%u') ...
      '_mrs' num2str(model.kv.mrsmsp,'%u') ...
      '_do' num2str(model.kv.do,'%u') ...
      '_tem' num2str(model.kv.SPLtem,'%u') 'dB_' model.flags.fbank];
if model.flags.do_gammatone
  cachename = [cachename '_'  num2str(1/model.kv.space,'%u') 'bpERB'];
  if model.flags.do_middleear; cachename = [cachename '_middleear']; end
  if model.flags.do_ihc; cachename = [cachename '_ihc']; end
else % zilany
  cachename = [cachename '_' model.flags.fibertypeseparation];
end
if model.kv.prior > 0 
  cachename = [cachename '_prior' num2str(model.kv.prior,'%u')];
end
if model.kv.tiwin < 0.5
  cachename = [cachename '_tiwin' num2str(model.kv.tiwin*1e3) 'ms']; 
end
cachename = [cachename '_mgs' num2str(model.kv.mgs)]; 


%% Hearing thresholds following OHC dysfunction
if flags.do_hearingthreshold

  cOHC = [1,0.5,0.1,0];
  flow = 125;%700; % Hz
  fhigh = 18000; % Hz
  spl = -20:100; % dB
  Ncf = 40; % # CF

  fs = 100e3; % Hz
  t = 0:1/fs:0.1; % s
  cf = audspace(flow,fhigh,Ncf); % CFs under test

  afr = amtcache('get',['hearingthreshold_' model.flags.fbank '_' model.flags.fibertypeseparation]);
  if isempty(afr)
    afr = nan(length(cf),length(cOHC),length(spl),3);
    for ff = 1:length(cf)
      sig = sin(2*pi*cf(ff)*t);
      for iiOHC = 1:length(cOHC)
        for iispl = 1:length(spl)
          for ft = 1:3
            ANoutTem = zilany2014(spl(iispl),sig,fs,...
                      'flow',cf(ff),'fhigh',cf(ff),...
                      'nfibers',1,'fiberType',ft,...
                      'cohc',cOHC(iiOHC),'cihc',1);
            afr(ff,iiOHC,iispl,ft) = mean(ANoutTem,2);
          end
        end
      end
      amtdisp([num2str(ff) ' of ' num2str(length(cf)) '  done.'],'progress')
    end
  amtcache('set',['hearingthreshold' cachename],afr)
  end


  ftd = [0.16,0.23,0.61]; % fiber type distribution from Liberman (1978)

  %% Evaluate avg. firing rate at normal hearing threshold (afrLMHc1)
  sizeafr = size(afr);
  afrLMH = sum(afr.*repmat(shiftdim(ftd,-2),[sizeafr(1:3),1]),4);
  afrLMHc1 = nan(length(cf),1);
  HT0 = nan(length(cf),1);
  for ff = 1:length(cf)
    [tmp,id] = min(abs(spl-absolutethreshold(cf(ff),'hda200')));
    HT0(ff) = spl(id);
    afrLMHc1(ff) = afrLMH(ff,1,id);
  end

  %% 
  % ftd(1:2) = 0;
  afr = sum(afr.*repmat(shiftdim(ftd,-2),[sizeafr(1:3),1]),4);
  HT = nan(length(cf),length(cOHC));
  afr0 = nan(length(cf),1);
  for ff = 1:length(cf)
    for iiOHC = 1:length(cOHC)
      id = find(afr(ff,iiOHC,:)>=afrLMHc1(ff),1,'first');
      if isempty(id)
        id = length(spl);
      end
      HT(ff,iiOHC) = spl(id);
    end
  end

  % Pure-tone average HL
  PTAf{1} = [500,1000,2000]; % Rakerd et al. (1998) 
  PTAf{2} = [3150,5000,8000]; % Rakerd et al. (1998) 
  PTAf{3} = [4000,8000,11000]; % Otte et al. (2013) 
  for cc = 1:length(cOHC)
    for ii = 1:length(PTAf)
      for ff = 1:length(PTAf{ii})
        HL(ff,cc,ii) = interp1(cf,HT(:,cc),PTAf{ii}(ff));
      end
    end
  end
  HL = HL(:,2:end,:) - repmat(HL(:,1,:),1,length(cOHC)-1,1);
  PTA = shiftdim(mean(HL));
  legendentries = cat(2,repmat('C_{OHC} = ',length(cOHC),1),num2str(cOHC(:),'%2.1f'));
  for cc = 1:length(cOHC)-1
    RowNames{cc} = legendentries(cc+1,:);
  end
  PTA = round(PTA);
  table(PTA(:,1),PTA(:,2),PTA(:,3),'RowNames',RowNames,'VariableNames',{'PTAlow_Rakerd98','PTAhigh_Rakerd98','PTAhigh_Otte13'})

  varargout{1} = HT;
  varargout{2} = cf;
  
  if flags.do_plot

    figure
    symb = {'-k*','-kd','-k>','-kp'};
    for ii = 1:size(HT,2)
      h(ii) = semilogx(cf,HT(:,ii),symb{ii});
      hold on
    end
    set(h,'MarkerFaceColor','w','MarkerSize',kv.MarkerSize)
    set(gca,'XLim',[cf(1),cf(end)],'YDir','reverse','FontSize',kv.FontSize,'TickLength',kv.TickLength)
    % set(gca,'XTickLabel',{'0.2','0.3','','0.5','','0.7','','','1','2','3','','5','','7','','','10'})
    set(gca,'Layer', 'top')
    xlabel('Frequency (Hz)','FontSize',kv.FontSize)
    ylabel('Hearing threshold (dB)','FontSize',kv.FontSize)

    leg = legend(legendentries,'Location','southwest');
    set(leg,'FontSize',kv.FontSize)

  end

end

%% ------ BASELINE EXAMPLES --------------------------------------------------------
if flags.do_baseline_ex
  
  latseg = 0;ii=1;%[-20,0,20]; ii = 2; % centers of lateral segments
%   dlat =  10;         % lateral range (+-) of each segment

  s = data_baumgartner2016('argimport',model.flags,model.kv);
  
%   idselect = ismember({s.id},{'NH15','NH22','NH62','NH12','NH39','NH18'});
  idselect = 19:23;
  s = s(idselect);

  %% LocaMo
  qe = zeros(length(s),length(latseg));
  pe = zeros(length(s),length(latseg));
  for ll = 1:length(s)

      s(ll).sphrtfs{ii} = 0;     % init
      s(ll).p{ii} = 0;        % init

      [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
      [s(ll).p{ii},respangs] = baumgartner2016(...
          s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},'argimport',model.flags,model.kv,...
          'ID',s(ll).id,'fs',s(ll).fs,...
          'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latseg(ii),'polsamp',polang,...
          'priordist',s(ll).priordist); 

      [ qe(ll,ii),pe(ll,ii) ] = baumgartner2014pmv2ppp( ...
          s(ll).p{ii} , polang , respangs , s(ll).target{ii});

      if flags.do_plot
        if ll ==1; figure; end
        subplot(2,3,ll)
        Nmax = min(150,s(ll).Ntar{ii});
        idplot = round(1:s(ll).Ntar{ii}/Nmax:s(ll).Ntar{ii});
        plot_baumgartner2014(s(ll).p{ii},polang,respangs,...
                  s(ll).target{ii}(idplot),s(ll).response{ii}(idplot),...
                  'MarkerSize',kv.MarkerSize,'cmax',0.05,'nocolorbar');
        title({['A: PE = ' num2str(s(ll).pe_exp_lat(ii),2) '\circ, QE = ' num2str(s(ll).qe_exp_lat(ii),2) '%'];['P: PE = ' num2str(pe(ll,ii),2) '\circ, QE = ' num2str(qe(ll,ii),2) '%']},...
          'FontSize',kv.FontSize-1)
        text(90,240,s(ll).id,'FontSize',kv.FontSize,...
          'Color','w','HorizontalAlignment','center')
        xlabel('Target Angle (deg)','FontSize',kv.FontSize)
        ylabel('Response Angle (deg)','FontSize',kv.FontSize)
        set(gca,'FontSize',kv.FontSize-1)
        set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
        set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})
      end

  end
  
  s = rmfield(s,{'Obj','itemlist','sphrtfs'}); % reduce file size 
  
  varargout{1} = s;
  
end

%% ------ PARAMETRIZATION -----------------------------------------------------------
if flags.do_parametrization
  
  if flags.do_mrsandgamma
  
    [gamma,mrs] = amtcache('get','parametrization', flags.cachemode);
    if isempty(gamma)
      amtdisp('Note that this procedure lasts at least 2 hours!','progress')

      tempfn = fullfile(amtbasepath,'experiments','exp_baumgartner2014_parametrization'); % temporary folder
      mkdir(tempfn)

      s = data_baumgartner2016('argimport',model.flags,model.kv);
      [gamma,mrs] = baumgartner2016parametrization(s,'SPLtem',kv.SPLtem,...
        flags.fbank,flags.fibertypeseparation,'mgs',kv.mgs);

      amtcache('set','parametrization',gamma,mrs);
    end
    varargout{1} = {gamma,mrs};
  
  else
    
    gamma = 7:2:30;%[3:2:7,8:13,15:2:30];
    
    dtot = nan(size(gamma));
    for g = 1:length(gamma)
      d = exp_baumgartner2016('argimport',model.flags,model.kv,'baseline','gamma',gamma(g));
      dtot(g) = d.total;
    end
    [d_opt,id_opt] = min(dtot);
    gamma_opt = gamma(id_opt);
    varargout{1} = gamma_opt;
    
    if flags.do_plot
      figure; 
      plot(gamma,dtot)
      xlabel('Gamma')
      ylabel('Prediction deviation')
    end
    
  end
  
end

%% ------ BASELINE ----------------------------------------------------------
if flags.do_baseline
  
  cachename = ['baseline_tar' num2str(model.kv.SPL,'%u') 'dB_' cachename];
  [r,d,s] = amtcache('get',cachename,flags.cachemode);
  
  if isempty(r)

    s = data_baumgartner2016('argimport',model.flags,model.kv);

%     qe_pred = nan(size(qe_exp));
%     pe_pred = qe_pred;
    for ii = 1:length(s)
      [p,rang,tang] = baumgartner2016(s(ii).Obj,s(ii).Obj,'argimport',model.flags,model.kv,...
            'ID',s(ii).id,'Condition','baseline','mrsmsp',s(ii).mrs,'S',s(ii).S,'priordist',s(ii).priordist);
      [s(ii).qe_pred,s(ii).pe_pred] = baumgartner2014pmv2ppp(p,tang,rang,s(ii).itemlist(:,6));
    end
    
    qe_exp = cat(1,s.qe_exp);
    pe_exp = cat(1,s.pe_exp);
    qe_pred = cat(1,s.qe_pred);
    pe_pred = cat(1,s.pe_pred);
    r_qe = corrcoef(qe_exp,qe_pred);
    r.qe = r_qe(2);
    r_pe = corrcoef(pe_exp,pe_pred);
    r.pe = r_pe(2);

    % prediction residues
    Ntargets = cat(1,s.Ntar); % # of targets
    Ntargets = cat(1,Ntargets{:});
    relfreq = Ntargets/sum(Ntargets(:));
    sd_pe = (pe_pred-pe_exp).^2; % squared differences
    d.pe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
    sd_qe = (qe_pred-qe_exp).^2;
    d.qe = sqrt(relfreq(:)' * sd_qe(:));

    s = rmfield(s,{'Obj'});
    
    amtcache('set',cachename,r,d,s)
  end
  
  d.total = (d.pe/90 + d.qe/100) /2;
  
  if nargout >0; varargout{1} = d;	end
  if nargout >1; varargout{2} = r;	end
  if nargout >2; varargout{3} = s;	end
  
  amtdisp(['Corr. QE: ' num2str(r.qe,'%2.2f') ', PE: ' num2str(r.pe,'%2.2f') ', QE+PE: ' num2str((r.qe+r.pe)/2,'%2.2f')])
    
  if flags.do_plot
    
    qe_exp = cat(1,s.qe_exp);
    pe_exp = cat(1,s.pe_exp);
    qe_pred = cat(1,s.qe_pred);
    pe_pred = cat(1,s.pe_pred);
    
    figure
    subplot(121)
    minqe = min([qe_exp(:);qe_pred(:)])-2;
    maxqe = max([qe_exp(:);qe_pred(:)])+2;
    limqe = [minqe,maxqe];
    plot(limqe,limqe,'k--')
    hold on
    plot(qe_exp,qe_pred,'ko')
    axis equal
    axis([minqe maxqe minqe maxqe])

    xlabel('Actual QE','FontSize',kv.FontSize)
    ylabel('Predicted QE','FontSize',kv.FontSize)
    title(['e_{QE} = ' num2str(d.qe,'%0.1f') '% , r_{QE} = ' num2str(r.qe,'%0.2f')],...
        'FontSize',kv.FontSize)

    subplot(122)
    minpe = min([pe_exp(:);pe_pred(:)])-2;
    maxpe = max([pe_exp(:);pe_pred(:)])+2;
    limpe = [minpe,maxpe];
    plot(limpe,limpe,'k--')
    hold on
    plot(pe_exp,pe_pred,'ko')
    axis equal
    axis([minpe maxpe minpe maxpe])
    xlabel('Actual PE','FontSize',kv.FontSize)
    ylabel('Predicted PE','FontSize',kv.FontSize)
    title(['e_{PE} = ' num2str(d.pe,'%0.1f') '\circ , r_{PE} = ' num2str(r.pe,'%0.2f')],...
        'FontSize',kv.FontSize)

  end
  
end


%% ------ SPATSTRAT EXAMPLES -----------------------------------------------------------
if flags.do_spatstrat_ex
  
  latdivision = 0;  % lateral angle
  dlat = 15;

  % Experimental Settings
  Conditions = {'BB','LP','W'};


  %% Computations
  s = data_baumgartner2016('argimport',model.flags,model.kv);  
  s = s(ismember({s.id},'NH58'));
  amtdisp(['Listener: ' s.id])
  chance = [];
  for C = 1:length(Conditions)

    Cond = Conditions{C};

    %% Data

    % Experimental data
    data = data_majdak2013(Cond);
    for ll = 1:length(s)
      if sum(ismember({data.id},s(ll).id)) % participant ?
        s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx;
        for ii = 1:length(latdivision)
          latresp = s(ll).itemlist(:,7);
          idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
          mm2 = s(ll).itemlist(idlat,:);
          chance = [chance;mm2];
          s(ll).target{ii} = mm2(:,6); % polar angle of target
          s(ll).response{ii} = mm2(:,8); % polar angle of response
        end
      end
    end


    for ll = 1:length(s)
        for ii = 1:length(latdivision)
            s(ll).spdtfs{ii} = 0;     % init
            s(ll).polang{ii} = 0;   % init
            [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(...
              latdivision(ii),s(ll).Obj);

            if C == 1       % Learn 
                s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};
            elseif C == 2   % Dummy
                temp=amtload('baumgartner2014','spatstrat_lpfilter.mat');
                s(ll).spdtfs_c{ii} = filter(temp.blp,temp.alp,s(ll).spdtfs{ii});
            elseif C == 3   % Warped
                s(ll).spdtfs_c{ii} = warp_hrtf(s(ll).spdtfs{ii},s(ll).fs);
            end

        end
    end


    %% Run Model

    for ll = 1:length(s)
      qe = zeros(1,length(latdivision));
      pe = zeros(1,length(latdivision));
      qe_t = zeros(1,length(latdivision));
      pe_t = zeros(1,length(latdivision));
      for ii = 1:length(latdivision)

        switch Cond
          case 'BB'
            clabel = 'baseline';
          case 'LP'
            clabel = 'lowpassed';
          case 'W'
            clabel = 'warped';
        end
        [s(ll).p{ii},rang] = baumgartner2016(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},...
              'argimport',model.flags,model.kv,...
              'ID',s(ll).id,'Condition',clabel,'fs',s(ll).fs,...
              'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);
        respangs{ii} = rang;

        [ qe(ii),pe(ii) ] = baumgartner2014pmv2ppp(s(ll).p{ii} , s(ll).polang{ii} , rang);

        if sum(ismember({data.id},s(ll).id)) % if participant then actual targets
          [ qe_t(ii),pe_t(ii) ] = baumgartner2014pmv2ppp( ...
              s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

        end

      end

      % Model results of pool
      s(ll).qe_pool(C,1) = mean(qe); 
      s(ll).pe_pool(C,1) = mean(pe);

      if sum(ismember({data.id},s(ll).id)) % participant ?
        % Actual experimental results
        s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
        s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
        s(ll).Nt(C,1) = size(s(ll).itemlist,1);
        % Model results of participants (actual target angles)
        if length(latdivision) == 3
          s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
              qe_t(2)*length(s(ll).target{2}) + ...
              qe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
              pe_t(2)*length(s(ll).target{2}) + ...
              pe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
        else 
          s(ll).qe_part(C,1) = mean(qe_t);
          s(ll).pe_part(C,1) = mean(pe_t);
        end

        if flags.do_plot

          if C == 1; figure; end
          subplot(1,3,C)
          ii = find(latdivision==0);
          responses = [];
          targets = [];
          for jj = ii
            responses = [responses;s(ll).response{jj}];
            targets = [targets;s(ll).target{jj}];
          end
          plot_baumgartner2014(s(ll).p{ii},s(ll).polang{ii},rang,...
                targets,responses,'MarkerSize',kv.MarkerSize,'cmax',0.05,'nocolorbar')
          text(90,240,Cond,...
            'FontSize',kv.FontSize,'Color','w','HorizontalAlignment','center')
          Nt = length(targets);
          tmp.m = [zeros(Nt,5) targets(:) zeros(Nt,1) responses(:)];
          tmp.qe = localizationerror(tmp.m,'querrMiddlebrooks');
          tmp.pe = localizationerror(tmp.m,'rmsPmedianlocal');
          title({['A: PE = ' num2str(tmp.pe,2) '\circ, QE = ' num2str(tmp.qe,2) '%'];...
            ['P: PE = ' num2str(pe_t(ii),2) '\circ, QE = ' num2str(qe_t(ii),2) '%']},'FontSize',kv.FontSize-1)
          xlabel('Target Angle (deg)','FontSize',kv.FontSize)
          ylabel('Response Angle (deg)','FontSize',kv.FontSize)
          set(gca,'FontSize',kv.FontSize-1)
          set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
          set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})

        end
        
      end

    end

  end
  
  varargout{1} = s;

end


%% ------ SPATSTRAT -----------------------------------------------------------
if flags.do_spatstrat
  
  cachename = ['spatstrat_tar' num2str(model.kv.SPL,'%u') 'dB_' cachename];
  [r,d,s,act,pred] = amtcache('get',cachename,flags.cachemode);
  
  if isempty(r)
    
    latdivision = 0;%[-20,0,20];            % lateral angle
    dlat = 30;%10;

    % Experimental Settings
    Conditions = {'BB','LP','W'};

    %% Computations
    s = data_baumgartner2016('argimport',model.flags,model.kv);
    for C = 1:length(Conditions)

      Cond = Conditions{C};

      %% Data

      % Experimental data
      data = data_majdak2013(Cond);
      for ll = 1:length(s)
        if sum(ismember({data.id},s(ll).id)) % if actual participant
          s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx; 
          for ii = 1:length(latdivision)
            latresp = s(ll).itemlist(:,7);
            idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
            mm2 = s(ll).itemlist(idlat,:);
            s(ll).target{ii} = mm2(:,6); % polar angle of target
            s(ll).response{ii} = mm2(:,8); % polar angle of response
          end
        end
      end


      for ll = 1:length(s)
          for ii = 1:length(latdivision)
              s(ll).spdtfs{ii} = 0;     % init
              s(ll).polang{ii} = 0;   % init
              [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(...
                latdivision(ii),s(ll).Obj);

              if C == 1       % Learn 
                  s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};
              elseif C == 2   % Dummy
                temp=amtload('baumgartner2014','spatstrat_lpfilter.mat');
                s(ll).spdtfs_c{ii} = filter(temp.blp,temp.alp,s(ll).spdtfs{ii});
              elseif C == 3   % Warped
                  s(ll).spdtfs_c{ii} = warp_hrtf(s(ll).spdtfs{ii},s(ll).fs);
              end

          end
      end


      %% Run Model
      idpart = [];
      for ll = 1:length(s)
        qe = zeros(1,length(latdivision));
        pe = zeros(1,length(latdivision));
        qe_t = zeros(1,length(latdivision));
        pe_t = zeros(1,length(latdivision));
        for ii = 1:length(latdivision)

          switch Cond
            case 'BB'
              clabel = 'baseline';
            case 'LP'
              clabel = 'lowpassed';
            case 'W'
              clabel = 'warped';
          end

          if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
            idpart = [idpart,ll];
            [s(ll).p{ii},rang] = baumgartner2016(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},'argimport',model.flags,model.kv,...
              'ID',s(ll).id,'Condition',clabel,'fs',s(ll).fs,...
              'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);
            
            [ qe_t(ii),pe_t(ii) ] = baumgartner2014pmv2ppp( ...
                s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

          end

        end

        if sum(ismember({data.id},s(ll).id)) % if actual participant 
          % Actual experimental results
          s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
          s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
          s(ll).Nt(C,1) = size(s(ll).itemlist,1);
          % Model results of participants (actual target angles)
          if length(latdivision) == 3
            s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
                qe_t(2)*length(s(ll).target{2}) + ...
                qe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
            s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
                pe_t(2)*length(s(ll).target{2}) + ...
                pe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          else 
            s(ll).qe_part(C,1) = mean(qe_t);
            s(ll).pe_part(C,1) = mean(pe_t);
          end

        end

      end
    end
    s = rmfield(s,{'Obj','spdtfs_c','spdtfs'});% reduce file size

    act.qe = [s(idpart).qe_exp]';
    act.pe = [s(idpart).pe_exp]';
    pred.qe = [s(idpart).qe_part]';
    pred.pe = [s(idpart).pe_part]';
    
    % Correlation coefficients
    [rmtx,p] =  corrcoef(act.qe,pred.qe);
    r.qe = rmtx(2);
    disp(['QE: r = ' num2str(r.qe,'%0.2f') ', p = ' num2str(p(2),'%0.3f')]);

    [rmtx,p] =  corrcoef(act.pe,pred.pe);
    r.pe = rmtx(2);
    disp(['PE: r = ' num2str(r.pe,'%0.2f') ', p = ' num2str(p(2),'%0.3f')]);


    % RMS Differences
    % individual:
    Ntargets = [s.Nt]'; % # of targets
    relfreq = Ntargets/sum(Ntargets(:));
    sd_pe = (pred.pe-act.pe).^2; % squared differences
    d.pe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
    sd_qe = (pred.qe-act.qe).^2;
    d.qe = sqrt(relfreq(:)' * sd_qe(:));
  
    amtcache('set',cachename,r,d,s,act,pred);
  end
  
  d.total = (d.pe/90 + d.qe/100) /2;
  
  if nargout >0; varargout{1} = d;	end
  if nargout >1; varargout{2} = r;	end
  if nargout >2; varargout{3} = s;	end
  
  if flags.do_plot
    
    %% Measures

    % Quartiles
    quart_pe_part = quantile(pred.pe,[.25 .50 .75]);
    quart_qe_part = quantile(pred.qe,[.25 .50 .75]);

    quart_pe_exp = quantile(act.pe,[.25 .50 .75]);
    quart_qe_exp = quantile(act.qe,[.25 .50 .75]);

    % Chance performance
    [qe0,pe0] = baumgartner2014pmv2ppp('chance');
    
    
    %% Plots
    dx = 0.15;
    
    figure 

    subplot(121)
    errorbar((1:3)+dx,quart_pe_part(2,:),...
        quart_pe_part(2,:) - quart_pe_part(1,:),...
        quart_pe_part(3,:) - quart_pe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');

    plot([0,4],[pe0,pe0],'k:')

    title(['e_{PE} = ' num2str(d.pe,'%0.1f') '\circ , r_{PE} = ' num2str(r.pe,'%0.2f')],...
      'FontSize',kv.FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',kv.FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[27 54.9],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    subplot(122)
    errorbar((1:3)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');
    l = legend('Model','Actual');
    set(l,'Location','northwest','FontSize',kv.FontSize-1)
    
    plot([0,4],[qe0 qe0],'k:')

    title(['e_{QE} = ' num2str(d.qe,'%0.1f') '% , r_{QE} = ' num2str(r.qe,'%0.2f')],...
      'FontSize',kv.FontSize)
    ylabel('Quadrant Error (%)','FontSize',kv.FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[0.1 54],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YAxisLocation','left',...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    set(gcf,'PaperPosition',[1,1,10,3.5])
    
  end
end

%% ------ NUMCHAN EXAMPLES -----------------------------------------------------------
if flags.do_numchan_ex
  
  % Model Settings
  latdivision = 0;            % lateral angle
  dlat = 10;

  % Experimental Settings
  Conditions = {'BB','CL','N24','N9','N3'};

  % Vocoder Settings 
  N = [inf,inf,24,9,3];
  flow = 300;     % lowest corner frequency
  fhigh = 16000;  % highest corner frequency

  %% Computations
  s = data_baumgartner2016('argimport',model.flags,model.kv);
  s = s(ismember({s.id},'NH33')); 
  disp(['Listener: ' s.id])
  chance = [];
  for C = 1:length(Conditions)

    Cond = Conditions{C};

    %% Data

    % Experimental data
    data = data_goupell2010(Cond);
    for ll = 1:length(s)
      if sum(ismember({data.id},s(ll).id)) % if actual participant
        s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx; 
        for ii = 1:length(latdivision)
          latresp = s(ll).itemlist(:,7);
          idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
          mm2 = s(ll).itemlist(idlat,:);
          chance = [chance;mm2];
          s(ll).target{ii} = mm2(:,6); % polar angle of target
          s(ll).response{ii} = mm2(:,8); % polar angle of response
        end
      end
    end

    % SP-DTFs
    for ll = 1:length(s)
        for ii = 1:length(latdivision)
            s(ll).spdtfs{ii} = 0;   % init
            s(ll).polang{ii} = 0;   % init
            [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(latdivision(ii),s(ll).Obj);
        end
    end


    %% Genereate conditional HRIRs

    stimPar.SamplingRate = s(ll).fs;
    imp = [1;zeros(2^12-1,1)]; % smooth results for 2^12
    for ll = 1:length(s)
      for ii = 1:length(latdivision)

        if strcmp(Cond,'BB') || strcmp(Cond,'CL')
            s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};

        else
          n = N(C);

          [syncrnfreq, GETtrain] = GETVocoder('',imp,n,flow,fhigh,0,100,stimPar);
          corners = [syncrnfreq(1);syncrnfreq(:,2)];

          ref = s(ll).spdtfs{ii};

          cond = zeros(length(imp),size(ref,2),2);

          for ch = 1:size(ref,3)
              for ang = 1:size(ref,2)
                  cond(:,ang,ch) = channelize('', 0.5*ref(:,ang,ch), ref(:,1), imp, n, corners, [], ...
                                  GETtrain, stimPar, 1, 0.01*s(ll).fs, 0.01*s(ll).fs);
              end
          end

          s(ll).spdtfs_c{ii} = cond;

        end
      end
    end

    if flags.do_plot 
     if C==1; fp = figure; end
     subplot(2,length(Conditions),C)
      if strcmp(Cond,'CL')
        stim = repmat([1;zeros(s(ll).fs/100-1,1)],17,1); % pulse train with 100pps
      else
        stim = noise(8e3,1,'white');
      end
      sig = lconv(stim,s(ll).spdtfs_c{ii});
     [mp,fc] = baumgartner2016spectralanalysis(sig,kv.SPL,...
       'argimport',model.flags,model.kv,'target','ID',s(ll).id,'Condition',Cond);
     pcolor(fc,s(ll).polang{ii},mp(:,:,1)')
     shading flat
     xlabel('Frequency (Hz)')
     ylabel('Discharge rate')
     title(Cond)
    end
    

    %% Run Model

    for ll = 1:length(s)
      clear qe pe qe_t pe_t
      for ii = 1:length(latdivision)
        if strcmp(Cond,'CL')
          stim = repmat([1;zeros(s(ll).fs/100-1,1)],10,1); % pulse train with 100pps
        else
          stim = [];
        end
        [p,rang] = baumgartner2016(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},'argimport',model.flags,model.kv,...
              'ID',s(ll).id,'Condition',Cond,...
              'stim',stim,'fsstim',s(ll).fs,...
              'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);

        [ qe(ii),pe(ii) ] = baumgartner2014pmv2ppp(p , s(ll).polang{ii} , rang);

        if sum(ismember({data.id},s(ll).id)) % if participant then actual targets
          [ qe_t(ii),pe_t(ii) ] = baumgartner2014pmv2ppp( ...
              p , s(ll).polang{ii} , rang , s(ll).target{ii} );
        end

      end

      % Model results of pool
      s(ll).qe_pool(C,1) = mean(qe); 
      s(ll).pe_pool(C,1) = mean(pe);

      if sum(ismember({data.id},s(ll).id)) % if actual participant
        % Actual experimental results
        s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
        s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
        s(ll).Nt(C,1) = size(s(ll).itemlist,1);
        % Model results of participants (actual target angles)
        if length(latdivision) == 3
          s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
              qe_t(2)*length(s(ll).target{2}) + ...
              qe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
              pe_t(2)*length(s(ll).target{2}) + ...
              pe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
        else 
          s(ll).qe_part(C,1) = mean(qe_t);
          s(ll).pe_part(C,1) = mean(pe_t);
        end

        if flags.do_plot && latdivision(ii) == 0

          subplot(2,length(Conditions),C+length(Conditions))

          plot_baumgartner2014(p,s(ll).polang{ii},rang,...
                s(ll).target{ii},s(ll).response{ii},...
                    'MarkerSize',kv.MarkerSize,'cmax',0.05,'nocolorbar');
          title({['A: PE = ' num2str(s(ll).pe_exp(C,1),2) '\circ, QE = ' num2str(s(ll).qe_exp(C,1),2) '%'];...
            ['P: PE = ' num2str(s(ll).pe_part(C,1),2) '\circ, QE = ' num2str(s(ll).qe_part(C,1),2) '%']},'FontSize',kv.FontSize-1)
          text(90,240,Cond,...
            'FontSize',kv.FontSize,'Color','w','HorizontalAlignment','center')
          xlabel('Target Angle (deg)','FontSize',kv.FontSize)
          ylabel('Response Angle (deg)','FontSize',kv.FontSize)
          set(gca,'FontSize',kv.FontSize-1)
          set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
          set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})

        end

      end

    end
  end
  
  varargout{1} = s;
  
end

%% ------ FIG NUMCHAN ----------------------------------------------------------
if flags.do_numchan
  
  cachename = ['numchan_tar' num2str(model.kv.SPL,'%u') 'dB_' cachename];
  [N,r,d,s] = amtcache('get',cachename,flags.cachemode);
  if isempty(N)
    
    % Model Settings
    latdivision = 0; % lateral angle
    dlat = 30;

    % Experimental Settings
    Conditions = {'BB','CL','N24','N18','N12','N9','N6','N3'};

    % Vocoder Settings 
    N = fliplr([3,6,9,12,18,24,30,36]);	% # of vocoder channels
    flow = 300;     % lowest corner frequency
    fhigh = 16000;  % highest corner frequency


    %% Computations
    s = data_baumgartner2016('argimport',model.flags,model.kv);
    for C = 1:length(Conditions)

      Cond = Conditions{C};

      %% Data

      % Experimental data
      data = data_goupell2010(Cond);
      for ll = 1:length(s)
        if sum(ismember({data.id},s(ll).id)) % if actual participant
          s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx; 
          for ii = 1:length(latdivision)
            latresp = s(ll).itemlist(:,7);
            idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
            mm2 = s(ll).itemlist(idlat,:); 
%             chance = [chance;mm2];
            s(ll).target{ii} = mm2(:,6); % polar angle of target
            s(ll).response{ii} = mm2(:,8); % polar angle of response
          end
        end
      end

      % SP-DTFs
      for ll = 1:length(s)
        for ii = 1:length(latdivision)
          s(ll).spdtfs{ii} = 0;   % init
          s(ll).polang{ii} = 0;   % init
          [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(latdivision(ii),s(ll).Obj);
        end
      end


      %% Genereate conditional HRIRs

      stimPar.SamplingRate = s(ll).fs;
      imp = [1;zeros(2^12-1,1)]; % smooth results for 2^12
      for ll = 1:length(s)
        for ii = 1:length(latdivision)

          if strcmp(Cond,'BB') || strcmp(Cond,'CL')
            s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};

          else
            n = N(C);
            cachenameGET = ['numchan_GET_N' num2str(n) '_' s(ll).id];
            cond = amtcache('get',cachenameGET);
            if isempty(cond)
              
              [syncrnfreq, GETtrain] = GETVocoder('',imp,n,flow,fhigh,0,100,stimPar);
              corners = [syncrnfreq(1);syncrnfreq(:,2)];

              ref = s(ll).spdtfs{ii};

              cond = zeros(length(imp),size(ref,2),2);

              for ch = 1:size(ref,3)
                for ang = 1:size(ref,2)
                    cond(:,ang,ch) = channelize('', 0.5*ref(:,ang,ch), ref(:,1), imp, n, corners, [], ...
                                    GETtrain, stimPar, 1, 0.01*s(ll).fs, 0.01*s(ll).fs);
                end
              end
              
              amtcache('set',cachenameGET,cond);
            end
            s(ll).spdtfs_c{ii} = cond;

          end
        end
      end

      %% Run Model

      idpart = [];
      for ll = 1:length(s)
        clear qe pe qe_t pe_t
        for ii = 1:length(latdivision)

          if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
            if strcmp(Cond,'CL')
              stim = repmat([1;zeros(s(ll).fs/100-1,1)],17,1); % pulse train with 100pps
            else
              stim = [];
            end
            
            [p,rang] = baumgartner2016(...
                s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},...
                'argimport',model.flags,model.kv,...
                'ID',s(ll).id,'Condition',Cond,...
                'stim',stim,'fsstim',s(ll).fs,...
                'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
                'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);
              
            [ qe_t(ii),pe_t(ii) ] = baumgartner2014pmv2ppp( ...
                p , s(ll).polang{ii} , rang , s(ll).target{ii} );
          end

        end

        if sum(ismember({data.id},s(ll).id)) % if actual participant 
          idpart = [idpart,ll];
          % Actual experimental results
          s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
          s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
          s(ll).Nt(C,1) = size(s(ll).itemlist,1);
          % Model results of participants (actual target angles)
          if length(latdivision) == 3
            s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
                qe_t(2)*length(s(ll).target{2}) + ...
                qe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
            s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
                pe_t(2)*length(s(ll).target{2}) + ...
                pe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          else 
            s(ll).qe_part(C,1) = mean(qe_t);
            s(ll).pe_part(C,1) = mean(pe_t);
          end
          
        end

      end
      disp(['Condition ' Cond ' completed.'])
    end
    
    act.qe = [s(idpart).qe_exp]';
    act.pe = [s(idpart).pe_exp]';
    pred.qe = [s(idpart).qe_part]';
    pred.pe = [s(idpart).pe_part]';
    
    [rmtx,p] =  corrcoef(act.qe,pred.qe);
    r.qe = rmtx(2);
%     r.qe.p = p(2);
    
    disp(['QE: r = ' num2str(r.qe,'%0.2f') ', p = ' num2str(p(2),'%0.3f')]);

    [rmtx,p] =  corrcoef(act.pe,pred.pe);
    r.pe = rmtx(2);
%     r.pe.p = p(2);
    disp(['PE: r = ' num2str(r.pe,'%0.2f') ', p = ' num2str(p(2),'%0.3f')]);
    
    % RMS Differences
    % individual:
    Ntargets = [s.Nt]'; % # of targets
    relfreq = Ntargets/sum(Ntargets(:));
    sd_pe = (pred.pe-act.pe).^2; % squared differences
    d.pe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
    sd_qe = (pred.qe-act.qe).^2;
    d.qe = sqrt(relfreq(:)' * sd_qe(:));
    
    s = rmfield(s,{'spdtfs','spdtfs_c','Obj','itemlist'});
    
    amtcache('set',cachename,N,r,d,s)
  end
  
  d.total = (d.pe/90 + d.qe/100) /2;
  
  if nargout >0; varargout{1} = d;	end
  if nargout >1; varargout{2} = r;	end
  if nargout >2; varargout{3} = s;	end
  if nargout >3; varargout{4} = N;	end
  

  if flags.do_plot
    
    data = data_goupell2010;
    idpart = ismember({s.id},{data.id});
    
    %% Quartiles
    quart_pe_part = fliplr(quantile([s(idpart).pe_part]',[.25 .50 .75]));
    quart_qe_part = fliplr(quantile([s(idpart).qe_part]',[.25 .50 .75]));

    quart_pe_exp = fliplr(quantile([s(idpart).pe_exp]',[.25 .50 .75]));
    quart_qe_exp = fliplr(quantile([s(idpart).qe_exp]',[.25 .50 .75]));

    [qe0,pe0] = baumgartner2014pmv2ppp('chance');
    
    
    %% Plot
    dx = 0.7;
    figure

    %% PE
    subplot(121)
    errorbar(fliplr(N)+dx,quart_pe_part(2,:),...
        quart_pe_part(2,:) - quart_pe_part(1,:),...
        quart_pe_part(3,:) - quart_pe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar(fliplr(N),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');
    plot([0,2*max(N)],[pe0,pe0],'k--')
    xlabel('Num. of channels','FontSize',kv.FontSize)
    ylabel('LPE (deg)','FontSize',kv.FontSize)

    title(['e = ' num2str(d.pe,'%0.1f') '\circ , r = ' num2str(r.pe,'%0.2f')],...
      'FontSize',kv.FontSize,'FontWeight','normal')
    set(gca,'XLim',[1 max(N)+2],'XTick',fliplr(N),...%[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;'';12;18;24;'CL';'BB'},...
        'YLim',[27 54.9],...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    %% QE
    subplot(122)
    errorbar(fliplr(N)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar(fliplr(N),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');
      
    l = legend('Model','Actual');
    set(l,'Location','northeast','FontSize',kv.FontSize-1)
      
    plot([0,2*max(N)],[qe0,qe0],'k--')
      
    title(['e = ' num2str(d.qe,'%0.1f') '% , r = ' num2str(r.qe,'%0.2f')],...
      'FontSize',kv.FontSize,'FontWeight','normal')
    xlabel('Num. of channels','FontSize',kv.FontSize)
    ylabel('QER (%)','FontSize',kv.FontSize)
    set(gca,'XLim',[1 max(N)+2],'XTick',fliplr(N),...%[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;'';12;18;24;'CL';'BB'},...
        'YLim',[0.1 54],...
        'YMinorTick','on',...
        'YAxisLocation','left','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
      

    set(gcf,'PaperPosition',[1,1,10,3.5])

  end
end

if flags.do_sabin2005
  
  cachename = ['sabin2005_' cachename '_SL2SPL' num2str(kv.SL2SPL) 'dB'];
  if model.kv.lat ~= 0
    cachename = [cachename '_lat' num2str(model.kv.lat)];
  end
  if flags.do_TLisSPL
    cachename = [cachename '_TLisSPL'];
  end
  if model.flags.do_gammatone
    cachename = [cachename '_minSPL' num2str(model.kv.GT_minSPL)];
  end
  pred = amtcache('get',cachename,flags.cachemode);
  if isempty(pred)
  
    s = data_baumgartner2016('argimport',model.flags,model.kv);
    
    SPL = [0:5:20,30:10:70]+kv.SL2SPL;

    pred.pvfront = nan(length(SPL),length(s));
    pred.pvrear = nan(length(SPL),length(s));
    pred.gfront = nan(length(SPL),length(s));
    pred.grear = nan(length(SPL),length(s));
    pred.precfront = nan(length(SPL),length(s));
    pred.precrear = nan(length(SPL),length(s));
    pred.prob = cell(length(SPL),length(s));
    for ii = 1:length(s)
% if ii == 1; figure; end
      for ll = 1:length(SPL)
        if flags.do_TLisSPL
          kv.SPLtem = SPL(ll);
        end
        if flags.do_nostat
          [pred.qe(ll,ii),pred.prob{ll,ii},m] = baumgartner2016(...
            s(ii).Obj,s(ii).Obj,...
            'argimport',model.flags,model.kv,...
            'ID',s(ii).id,'S',s(ii).S,'SPL',SPL(ll),...
            'priordist',s(ii).priordist,'QE');
          pred.pvfront(ll,ii) = nan;
          pred.pvrear(ll,ii) = nan;
          pred.gfront(ll,ii) = nan;
          pred.grear(ll,ii) = nan;
          pred.precfront(ll,ii) = nan;
          pred.precrear(ll,ii) = nan;
        else
          [pred.pvfront(ll,ii),pred.prob{ll,ii},m] = baumgartner2016(...
            s(ii).Obj,s(ii).Obj,...
            'argimport',model.flags,model.kv,...
            'ID',s(ii).id,'S',s(ii).S,'SPL',SPL(ll),...
            'priordist',s(ii).priordist,'pVeridicalPfront');
          pred.pvrear(ll,ii) = localizationerror(m,'pVeridicalPrear');
          pred.gfront(ll,ii) = localizationerror(m,'gainPfront');
          pred.grear(ll,ii) = localizationerror(m,'gainPrear');
          pred.precfront(ll,ii) = localizationerror(m,'precPregressFront');
          pred.precrear(ll,ii) = localizationerror(m,'precPregressRear');
        end
% if ii == 1
%   subplot(3,3,ll)
%   plot_baumgartner2014(pred.prob{ll,ii}.p,pred.prob{ll,ii}.tang,pred.prob{ll,ii}.rang,m(:,6),m(:,8))
%   colorbar off
%   title(SPL(ll))
% end
      end
      amtdisp([num2str(ii,'%u') ' of ' num2str(length(s),'%u') ' done'],'progress')
    end

    pred.SPL = SPL;
    amtcache('set',cachename,pred)
  end
  
  data = data_sabin2005;
  data.SPL = data.SL + kv.SL2SPL; % assumption on SL
  
  %% Gain conversion: slope in deg -> limited range and equidistant
  if flags.do_gainInDeg
    pred.gfront = gain2slope(pred.gfront);
    pred.grear = gain2slope(pred.grear);
    data.gain.f.m = gain2slope(data.gain.f.m);
    data.gain.r.m = gain2slope(data.gain.r.m);
    data.gain.f.sd = gain2slope(data.gain.f.sd);
    data.gain.r.sd = gain2slope(data.gain.r.sd);
  end
  
  %% Restriction to reliable data (at least 75% audible trials in Sabin2005)
  dvar = {'gain','pqv','var'}; % data variable names
  mvar = {'g','pv','prec'}; % modeled variable names
  SL = pred.SPL-kv.SL2SPL; % assumption on SL
  SPL = pred.SPL;
  
  minSLf = 10;
  minSLr = 15;
  
  for ii = 1:length(mvar)
    eval(['data.' dvar{ii} '.f.m(data.SL < minSLf) = nan;'])
    eval(['data.' dvar{ii} '.f.sd(data.SL < minSLf) = nan;'])
    eval(['pred.' mvar{ii} 'front(SL < minSLf,:) = nan;'])
    
    eval(['data.' dvar{ii} '.r.m(data.SL < minSLr) = nan;'])
    eval(['data.' dvar{ii} '.r.sd(data.SL < minSLr) = nan;'])
    eval(['pred.' mvar{ii} 'rear(SL < minSLr,:) = nan;'])
  end
  
  %% Prediction deviation score
  limits = {'90','100','45'}; % for normalization
  idc = ismember(pred.SPL,data.SPL);
  
  d = zeros(length(mvar),2);
  for ii = 1:length(mvar)
    eval(['d20pooled.f = [data.' dvar{ii} '.f.m(1:4),mean(data.' dvar{ii} '.f.m(5:6)),data.' dvar{ii} '.f.m(7:end)];']) % 20dB SL pooled
    eval(['d20pooled.r = [data.' dvar{ii} '.r.m(1:4),mean(data.' dvar{ii} '.r.m(5:6)),data.' dvar{ii} '.r.m(7:end)];'])
    if flags.do_nostat
      eval(['d(ii,1) = sqrt(mean((mean(pred.' mvar{ii} 'front(idc,:)),2) - transpose(d20pooled.f)).^2))/' limits{ii} ';'])
      eval(['d(ii,2) = sqrt(mean((mean(pred.' mvar{ii} 'rear(idc,:),2) - transpose(d20pooled.r)).^2))/' limits{ii} ';'])
    else
      eval(['d(ii,1) = sqrt(nanmean((nanmean(pred.' mvar{ii} 'front(idc,:),2) - transpose(d20pooled.f)).^2))/' limits{ii} ';'])
      eval(['d(ii,2) = sqrt(nanmean((nanmean(pred.' mvar{ii} 'rear(idc,:),2) - transpose(d20pooled.r)).^2))/' limits{ii} ';'])
    end
  end
  d = mean(d(:));
  amtdisp(['Prediction deviation score: ' num2str(d)])
  
  %% Output
  varargout{1}=d;
  varargout{2}=pred;
  varargout{3}=data;

  %% Plot
  if flags.do_plot
    
    minSPLf = minSLf + kv.SL2SPL;
    minSPLr = minSLr + kv.SL2SPL;
    maxSPL = 70 + kv.SL2SPL;
    marSPL = 4.9; % margin of SL
    
    if flags.do_gainInDeg
      ylim = {[-2,60],[-5,105],[8,42]}; % ylimits acc. to Sabin et al. (2005)
      elabel = {'Slope (deg)','% Quasi-veridical','Variability (deg)'}; % plot ylabels
    else
      ylim = {[-0.2,1.6],[-5,105],[8,42]}; % ylimits acc. to Sabin et al. (2005)
      elabel = {'Gain','% Quasi-veridical','Variability (deg)'}; % plot ylabels
    end
    
    figure;
    ha = tight_subplot(3,2,kv.gap,kv.marg_h,kv.marg_w);
    for ii = 1:length(mvar)
%       subplot(3,2,1+(ii-1)*2)
      axes(ha(1+(ii-1)*2))
      if ii == 1
        plot([-10,100],[45,45],'k--') % ideal slope
        hold on
      end
      eval(['p1 = errorbar(SPL-0.3,nanmean(pred.' mvar{ii} 'front,2),nanstd(pred.' mvar{ii} 'front,1,2));'])
      hold on
      eval(['p2 = errorbar(data.SPL+0.3,data.' dvar{ii} '.f.m,data.' dvar{ii} '.f.sd);'])
      set(p1,'MarkerFaceColor','k','Marker','^','Color','k','MarkerSize',kv.MarkerSize)
      set(p2,'MarkerFaceColor','w','Marker','^','Color','k','MarkerSize',kv.MarkerSize)
%       axis([-8,68,ylim{ii}])
      axis([minSPLf-marSPL,maxSPL+marSPL,ylim{ii}])
      ylabel(elabel{ii},'FontSize',kv.FontSize)
      set(gca,'XTick',round(minSPLf/10)*10:10:maxSPL)
      set(gca,'TickLength',2*get(gca,'TickLength'),'FontSize',kv.FontSize)

      if ii == 1
        title('Front','FontSize',kv.FontSize)
      elseif ii == 3
        xlabel('SPL (dB)','FontSize',kv.FontSize)
%         leg = legend('Model','Actual');
%         set(leg,'Location','north','FontSize',kv.FontSize)
      end
      if ii<3
        set(gca,'XTickLabel',[])
      end
      
%       subplot(3,2,2+(ii-1)*2)
      axes(ha(2+(ii-1)*2))  
      if ii == 1
        plot([-10,100],[45,45],'k--') % ideal slope
        hold on
      end    
      eval(['p1 = errorbar(SPL-0.3,nanmean(pred.' mvar{ii} 'rear,2),nanstd(pred.' mvar{ii} 'rear,1,2));'])
      hold on
      eval(['p2 = errorbar(data.SPL+0.3,data.' dvar{ii} '.r.m,data.' dvar{ii} '.r.sd);'])
      set(p1,'MarkerFaceColor','k','Marker','o','Color','k','MarkerSize',kv.MarkerSize)
      set(p2,'MarkerFaceColor','w','Marker','o','Color','k','MarkerSize',kv.MarkerSize)
%       axis([-8,68,ylim{ii}])
      axis([minSPLr-marSPL,maxSPL+marSPL,ylim{ii}])
      set(gca,'XTick',round(minSPLr/10)*10:10:maxSPL,'YTickLabel',[])
      set(gca,'TickLength',2*get(gca,'TickLength'),'FontSize',kv.FontSize)
%       ylabel(elabel{ii},'FontSize',kv.FontSize)
      if ii == 1
        title('Rear','FontSize',kv.FontSize)
      elseif ii == 3
        xlabel('SPL (dB)','FontSize',kv.FontSize)
%         leg = legend('Model','Actual');
%         set(leg,'Location','north','FontSize',kv.FontSize)
      end
      
    end
  end
  
end

if flags.do_impairment
  
  if isempty(errorflag)
    errorflag = 'QE';
    amtdisp('Localization performance measure not chosen -> QE used.')
  end
  
  cohc = [1,0.5,0.1,0];
  ft = {1:3;1;2;3};
  SPL = [80,50];
  
  cachename = ['impairment_' cachename '_' errorflag];
  s = amtcache('get',cachename,flags.cachemode);
  if isempty(s)
    
    s = data_baumgartner2016('argimport',model.flags,model.kv);
      
    if not(isempty(kv.subjects))
      s = s(kv.subjects);
    end

    for ii=1:length(s)
      for cc=1:length(cohc)
        for ll=1:length(SPL)
          for ff=1:length(ft)
            err = baumgartner2016(s(ii).Obj,s(ii).Obj,...
              'argimport',model.flags,model.kv,'ID',s(ii).id,errorflag,'fiberTypes',ft{ff},...
              'S',s(ii).S,'cohc',cohc(cc),'SPL',SPL(ll),'priordist',s(ii).priordist);
            s(ii).err(cc,ff,ll) = err;
            s(ii).cohc(cc,ff,ll) = cohc(cc);
            s(ii).SPL(cc,ff,ll) = SPL(ll);
            s(ii).ft{cc,ff,ll} = ft{ff};
          end
        end
      end
      amtdisp([num2str(ii,'%u') ' of ' num2str(length(s),'%u') ' done'],'progress')
    end
    
    s = rmfield(s,{'Obj','itemlist'});
    amtcache('set',cachename,s);
  end
  varargout{1} = s;
  
  err = cat(4,s.err);
  
  Ncohc = size(err,1);
  Nft = size(err,2);
  Nspl = size(err,3);
  Ncond = length(s(1).cohc(:)); % #conditions
  Nsub = length(s); % #subjects
  mtx.err = reshape(shiftdim(err,3),[Nsub,Ncond]);
  
  % labels
  cohcstr = num2str(s(1).cohc(:),'%2.1f');
  SPLstr = num2str(s(1).SPL(:),'%2.0f');
  ftnum = s(1).ft(:);
  ftstr = cell(length(ftnum),1);
  XTickLabel = ftstr;
  for ii = 1:length(ftnum)
    lab = ['cohc: ' cohcstr(ii,:) ', ' SPLstr(ii,:) 'dB, '];
    if ftnum{ii} == 1
      ftstr{ii} = 'LSR';
    elseif ftnum{ii} == 2
      ftstr{ii} = 'MSR';
    elseif ftnum{ii} == 3
      ftstr{ii} = 'HSR';
    else % ft{ii} == 1:3
      ftstr{ii} = 'All';
    end
    XTickLabel{ii} = [lab ftstr{ii}];
  end
    
  % sort data acc. to ascending SPL
  [tmp,idSPLsort] = sort(s(1).SPL(:));
  ftstr = ftstr(idSPLsort);
  cohcstr = cellstr(cohcstr(idSPLsort,:));
  SPLstr = cellstr(SPLstr(idSPLsort,:));
  mtx.err = mtx.err(:,idSPLsort,:);
    
  if flags.do_plot
    
    % interaction plots
    merr = mean(err,4); % average across subjects; dims: [cohc,ft,SPL]
    figure
    % OHC-SPL
    subplot(1,3,1)
    err_OHC_SPL = squeeze(mean(merr,2));
    plot(cohc,err_OHC_SPL)
    legend(num2str(SPL(:)))
    xlabel('C_{OHC}')
    ylabel(errorflag)
    % OHC-FT
    subplot(1,3,2)
    err_OHC_FT = squeeze(mean(merr,3));
    plot(cohc,err_OHC_FT)
    legend('All','LSR','MSR','HSR')
    xlabel('C_{OHC}')
    % SPL-FT
    subplot(1,3,3)
    err_FT_SPL = squeeze(mean(merr,1));
    plot(SPL,err_FT_SPL)
    legend('All','LSR','MSR','HSR')
    xlabel('SPL')
    
    
    emax = max(mtx.err(:))+0.5;
    emin = min(mtx.err(:))-0.5;
    De = emax-emin;
    
    figure;
    b = boxplot(mtx.err,{ftstr(:),cohcstr(:),SPLstr(:)},...
        'plotstyle','compact','colors',repmat([.5,.5,.5;0,0,0],Ncond/Nspl,1),'medianstyle','line',...
        'factorgap',[],'labelverbosity','all','symbol','');
    hold on
    if flags.do_FTlabel
      for ii = 1:Nft
        jj = 1+(ii-1)*Nspl*Ncohc;
        yftlbl = emax;
        text(jj,yftlbl,ftstr{ii*Ncohc},'FontSize',kv.FontSize) % inside panel
      end
    end
    
    [qe0,pe0] = baumgartner2014pmv2ppp('chance');
    if strcmp(errorflag,'QE')
      h = plot([0.5,Ncond+0.5],[qe0,qe0],'k--');uistack(h, 'bottom')
    elseif strcmp(errorflag,'PE')
      h = plot([0.5,Ncond+0.5],[pe0,pe0],'k--');uistack(h, 'bottom')
    end
    
    axis([0.5,Ncond+0.5,emin-De/20,emax+De/8])
    set(gca,'XTick',1.5:2:Ncond,'XTickLabel',cohcstr(1:Ncohc),'FontSize',kv.FontSize);
    xlabel({'C_{OHC}'},'FontSize',kv.FontSize,'FontWeight','bold')
    ylabel(errorflag,'FontSize',kv.FontSize,'FontWeight','bold') 
    
  end
  
  if flags.do_stat
    s = rmfield(s,{'pe_exp','qe_exp','S','pe_exp_lat','qe_exp_lat','target','response','Ntar','mrs','fs'});
    
    t = array2table(mtx.err);
    within = table(SPLstr,cohcstr,ftstr,'VariableNames',{'SPL','Cohc','FT'});
    rm = fitrm(t,['Var1-Var' num2str(Ncond) ' ~ 1'],'WithinDesign',within); % no between-subjects factors -> only intercept
    
    % 3-way repeated-measures ANOVA
    [tbl.ranova,A,C,D] = ranova(rm,'WithinModel','Cohc*SPL*FT');
    tbl.ranova.Properties.RowNames = strrep(tbl.ranova.Properties.RowNames,'(Intercept):','');
    
    % Mauchly's test for sphericity
    tbl.mauchly = mauchly(rm,C);
    tbl.mauchly.Properties.RowNames = tbl.ranova.Properties.RowNames(1:2:end);
    
    % Sphericity corrections
    tbl.eps = epsilon(rm,C);
    tbl.eps.Properties.RowNames = tbl.ranova.Properties.RowNames(1:2:end);
    
    % Add corrected DFs to ranova table
    idrep = round(0.5:0.5:length(tbl.eps.GreenhouseGeisser)); % repeat iteratively
    tbl.ranova.DFGG = tbl.ranova.DF .* tbl.eps.GreenhouseGeisser(idrep);
    
    % Post-hoc analysis
    tbl.posthoc.Cohc = multcompare(rm,'Cohc');
    tbl.posthoc.FT = multcompare(rm,'FT');
    
    % Display results
    amtdisp(['3-way repeated-measures ANOVA for ' errorflag])
    amtdisp(tbl.ranova)
    amtdisp('Mauchly test and sphericity corrections')
    amtdisp([tbl.mauchly,tbl.eps])
    amtdisp('Posthoc analysis')
    amtdisp(tbl.posthoc.Cohc)
    amtdisp(tbl.posthoc.FT)
    
    varargout{1} = tbl;
    
  end
end

if flags.do_effectOnCues
  
  sid = 10;    % listener No.
  spl = kv.SPL;   % SPL in dB
%   tang = 0;   % target polar angle
  
  s = data_baumgartner2016('argimport',model.flags,model.kv);
  [dtf,polang] = extractsp(0,s(sid).Obj);
  sig = lconv(noise(8e3,1,'white'),dtf);
      
  amtdisp(['Exemplary listener: ' s(sid).id])
  
  ymin = 0;
  ymax = kv.mgs*pi;
  spl = [50,80];
  
  if flags.do_FT
    % Effect of FT  
    FT = {1:3,1,2,3};
    mp_ft = cell(length(FT),length(spl));
    gp_ft = cell(length(FT),length(spl));
    for ll = 1:length(spl)
      for ft = 1:length(FT)
        [mp_ft{ft,ll},fc] = baumgartner2016spectralanalysis(sig,spl(ll),...
          'target','ID',s(sid).id ,'Condition','baseline','fiberTypes',FT{ft},flags.amtcache);
        [gp_ft{ft,ll},gfc] = baumgartner2016gradientextraction(mp_ft{ft,ll},fc,'mgs',kv.mgs);
      end
    end

    if flags.do_plot
      figure
      ha = tight_subplot(length(FT),length(cohc),kv.gap,kv.marg_h,kv.marg_w);
      colormap gray
      colormap(flipud(colormap))
      labels = {'LMH','L','M','H'};
      for ll = 1:length(spl)
        for ft = 1:length(FT)
          axes(ha(ft+(ll-1)*length(FT)))
          pcolor(gfc,polang,gp_ft{ft,ll}.m(:,:,1)')
          shading flat
          caxis([ymin,ymax])
          title(labels{ft})
          xlabel('Frequency (kHz)')
          ylabel('Polar angle (deg)')
          set(gca,'XScale','log','FontSize',kv.FontSize)
          set(gca,'layer','top',...
            'XTick',[1,2,4,8,16]*1e3,...
            'XTickLabel',[1,2,4,8,16],...
            'YTick',-30:30:180,...
            'FontSize',kv.FontSize)
        end
      end
    end

  else
    % Effect of C_OHC
    cohc = [1,0.5,0.1,0];
    gp_cohc = cell(length(cohc),length(spl));
    for ll = 1:length(spl)
      for cc = 1:length(cohc)
        [mp,fc] = baumgartner2016spectralanalysis(sig,spl(ll),'target',...
          'ID',s(sid).id,'Condition','baseline','cohc',cohc(cc),flags.amtcache);
        [gp_cohc{cc,ll},gfc] = baumgartner2016gradientextraction(mp,fc,'mgs',kv.mgs);
      end
    end

    if flags.do_plot
      figure
      ha = tight_subplot(length(spl),length(cohc),[.02,.01],kv.marg_h,kv.marg_w);
      colormap gray
      colormap(flipud(colormap))
      labels = {'C_{OHC} = 1','C_{OHC} = 0.5','C_{OHC} = 0.1','C_{OHC} = 0'};
      for ll = 1:length(spl)
        for cc = 1:length(cohc)
          axes(ha(cc+(ll-1)*length(cohc)))
          pcolor(gfc,polang,gp_cohc{cc,ll}.m(:,:,1)')
          shading flat
          caxis([ymin,ymax])
          % xlabel and COHC
          set(gca,'XTick',[1,2,4,8,16]*1e3,'XTickLabel',[1,2,4,8,16])
          if ll == 2
            if cc == 2
              xlabel(['                                            ',...
                'Frequency (kHz)'],'FontSize',kv.FontSize)
            end
          else
            title(labels{cc},'FontSize',kv.FontSize)
            set(gca,'XTickLabel',[])
          end
          % ylabel and SPL
          if cc == 1
            ylabel('Polar angle (deg)','FontSize',kv.FontSize)
            text(180,200,[num2str(spl(ll)) ' dB'],'FontWeight','bold','FontSize',kv.FontSize)
          else
            set(gca,'YTickLabel',[])
          end
          set(gca,'XScale','log','FontSize',kv.FontSize)
          set(gca,'layer','top',...
            'TickLength',2*get(gca,'TickLength'),... 
            'YTick',-30:30:180,...
            'FontSize',kv.FontSize)
        end
      end
      
      c = colorbar;
%       pos = get(c,'Position');
      set(c,'Position',[.935,.25,.02,.5])
      set(get(c,'Label'),'String','Spikes/s/ERB','FontSize',kv.FontSize)
      
    end
  end
  
end

if flags.do_dynrangecheck
  
  splmax = 140; % dB
  splHRTF = [50,80]; % dB

  panels = {'C_{OHC} = 1','C_{OHC} = 0.5','C_{OHC} = 0.1','C_{OHC} = 0'};
  if flags.do_separate
    labels = {'LSR','MSR','HSR'};
  else
    labels = {'LMH','MH','H'};
  end

  splminplot = 15;
  
  figure
  ha = tight_subplot(4,1,kv.gap,kv.marg_h,kv.marg_w);
  
  cohc = [1,0.5,0.1,0];
  for cc = 1:length(cohc)
  
    sig = noise(0.1*48e3,1,'white'); % 100-ms Gaussian white noise burst
    spl = 0:10:splmax;
    mp = zeros(kv.nf,length(spl),2,3);
    if flags.do_separate
      for ii = 1:length(spl)
        [mp(:,ii,:,:),fc] = baumgartner2016spectralanalysis(cat(3,sig,sig),spl(ii),...
            'target','ID','','Condition','dynrangecheck','fiberTypes',1:3,'ftopt','cohc',cohc(cc));
      end
    else
      fiberTypes = {1:3,2:3,3};
      for ft = 1:length(fiberTypes)
        for ii = 1:length(spl)
          [mp(:,ii,:,ft),fc] = baumgartner2016spectralanalysis(cat(3,sig,sig),spl(ii),...
              'target','dynrangecheck','fiberTypes',fiberTypes{ft},'cohc',cohc(cc));
        end
      end
    end

    if flags.do_dynrangeDiff
      mp = diff(mp,1,2)/mean(diff(spl));
      spl = spl(2:end);
    end
    
    if flags.do_dprime
      sd = 2.6*mp.^0.34;
      sdDiff = sqrt(sd(:,1:length(spl)-1,:,:).^2 + sd(:,2:length(spl),:,:).^2);
      mp = diff(mp,1,2)./sdDiff;
      spl = spl(2:end);
    end

    axes(ha(cc))
    
    % target HRTF range
    color = {.9*ones(1,3),.8*ones(1,3)};
    for ll = 1:length(splHRTF)
      a(1) = area(splHRTF(ll)+[-10,10],[1e3,1e3],'EdgeColor',ones(1,3));
      hold on
      a(2) = area(splHRTF(ll)+[-10,10],-[1e3,1e3],'EdgeColor',ones(1,3));
      set(a,'FaceColor',color{ll})
    end
    set(gca,'Layer','top')
    
    % Rate-intensity curves
    sty = {'kv-','ks-','k^-'};
    for ft = 1:3
      h(ft) = plot(spl,mean(mp(:,:,1,ft),1),sty{ft});
      hold on
      if cc == length(cohc)
        xlabel('SPL (dB)','FontSize',kv.FontSize)
      end
      if flags.do_dynrangeDiff
        text(splminplot+5,9,panels{cc},'FontSize',kv.FontSize)
        ylabel({'Rate difference (spikes/s/dB)'},'FontSize',kv.FontSize)
      elseif flags.do_dprime
        text(splminplot+5,5,panels{cc},'FontSize',kv.FontSize)
        if cc==3
          text(0,6,{'d^{\prime}'},'FontSize',kv.FontSize)
        end
%         ylabel({'d^{\prime}'},'FontSize',kv.FontSize)
      else
        text(splminplot+5,330,panels{cc},'FontSize',kv.FontSize)
        ylabel('Firing rate (spikes/s)','FontSize',kv.FontSize)
      end
    end
    set(h,'MarkerFaceColor','k')
    if cc == 1%length(cohc)
      leg = legend(h,labels);
      set(leg,'Location','northoutside','FontSize',kv.FontSize,'Orientation','horizontal')
      set(leg,'Position',[.4,.96,.33,.03])
    end
    
    if flags.do_dynrangeDiff
      plot([0,splmax],[0,0],'k--') % sensitivity threshold
      axis([splminplot,splmax-5,-1,10.5])
    elseif flags.do_dprime
      axis([splminplot,splmax-5,-1,5.9])  
    else
      axis([splminplot,splmax-5,-20,399])
    end
    XTick = 20:10:130;
    XTickLabel = num2cell(XTick);
    XTickLabel(2:2:end) = {' '};
    set(gca,'XTick',XTick,'XTickLabel',XTickLabel,'FontSize',kv.FontSize)
  
  end
  
end
  
if flags.do_localevel
  
  dlat = 30; % deg
  condition = data_baumgartner2015('ConditionNames');
  
  cachename = ['localevel_' cachename];
  if not(kv.gammashortfact == 1)
    cachename = [cachename '_gsf' num2str(kv.gammashortfact,'%1.1f')];
  end
  if not(kv.Sshortfact == 1)
    cachename = [cachename '_Ssf' num2str(kv.Sshortfact,'%1.1f')];
  end
  if not(kv.psgeshort == 1)
    cachename = [cachename '_psgec' num2str(kv.psgeshort,'%1.1f')];
  end
  if flags.do_TLisSPL
    cachename = [cachename '_TLisSPL'];
  end
  if flags.do_gammatone
    cachename = [cachename '_minSPL' num2str(kv.GT_minSPL) '_maxSPL' num2str(kv.GT_maxSPL)];
  end
  [pred,ref] = amtcache('get', cachename, flags.cachemode);
  if isempty(pred)
    pred.qe = nan(7,length(condition));
    pred.pe = pred.qe;
    ref.qe = pred.qe;
    ref.pe = pred.qe;
    for cc = 1:length(condition)

      data = data_baumgartner2015(condition{cc},'model','gamma',kv.gamma,...
        'tiwin',kv.tiwin,flags.fbank,flags.fibertypeseparation,flags.recalib,'mgs',kv.mgs);

      for ll = 1:length(data)

        latresp = data(ll).itemlist(:,7);
        idlat = latresp <= dlat & latresp > -dlat;
        itemlist = data(ll).itemlist(idlat,:);
        exptang = itemlist(:,6);
        exprang = itemlist(:,8);

        ref.pe(ll,cc) = real(localizationerror(itemlist,'rmsPmedianlocal'));
        ref.qe(ll,cc) = real(localizationerror(itemlist,'querrMiddlebrooks'));

        if flags.do_TLisSPL
          kv.SPLtem = data(ll).SPL;
        end
        if strcmp(condition{cc},'Long')
          clbl = condition{cc};
        else % short
          clbl = '';
        end
        [err,pmv] = baumgartner2016(data(ll).Obj,data(ll).Obj,...
            'argimport',model.flags,model.kv,'ID',data(ll).id,'Condition',clbl,'S',data(ll).S,...
            'stim',data(ll).stim,'fsstim',data(ll).fsstim,'SPL',data(ll).SPL,...
            'QE_PE_EB','exptang',exptang,'priordist',data(ll).priordist);
        pred.pmv{ll,cc} = pmv;
        pred.qe(ll,cc) = err.qe;
        pred.pe(ll,cc) = err.pe;

      end
    end
    amtcache('set',cachename,pred,ref);
  else
    data = data_baumgartner2016;
  end
  
  perrmtx = ref.pe';
  pred.perrmtx = pred.pe';
  qerrmtx = ref.qe';
  pred.qerrmtx = pred.qe';
  
  Nsub = length(data);
  
  %% Prediction residues 
  mm = 1;
  [r,p] = corrcoef([pred.pe],[ref.pe]);
  r_perr(mm) = r(2);
  e_perr(mm) = mean(rms([pred.pe]-[ref.pe]));

  [r,p] = corrcoef([pred.qe],[ref.qe]);
  r_qerr(mm) = r(2);
  e_qerr(mm) = mean(rms([pred.qe]-[ref.qe]));
    
  amtdisp(' e_PE  r_PE  e_QE   r_QE')
  amtdisp([num2str(e_perr(mm),'%2.1f') '\deg  ' num2str(r_perr(mm),'%2.2f') '  ' num2str(e_qerr(mm),'%2.1f') '%  ' num2str(r_qerr(mm),'%2.2f')])
  
  varargout{1} = 0.5* (e_perr(mm)/90 + e_qerr(mm)/100);
  
  %% Plots
  if flags.do_plot
    
    if flags.do_performance
    
      LineWidth = 1;

      symb = {...
            'rs-';...
            'bd-';...
            'gh-';...
            'mv-';...
            'c*-';...
            };
      symbExp = 'ko--';

      name{mm} = 'Pred.';
      % Legend
      legendentry = {'Actual, long';'Actual, short'};
    %   for mm = 1:length(pred)
        legendentry{2+2*mm-1} = [name{mm} ', long'];
        legendentry{2+2*mm} = [name{mm} ', short'];
    %   end

      xval = 10:10:70;
      xtext = 10;

      % Local central RMS error
      hfig = figure;
      ha = tight_subplot(8,2,kv.gap,kv.marg_h,kv.marg_w);
      for ii=1:Nsub

        axes(ha(2*ii));
        hlong = plot(50,perrmtx(1,ii),symbExp(1:2));
        set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)

        set(gca,'YAxisLocation','right','YMinorTick','on')
    %     set(gca,'XTick',xval,'YTick',10:20:50)
        set(gca,'TickLength',kv.TickLength)

        hold on
        hshort = plot(xval,perrmtx(2:end,ii),symbExp);
        set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)

        for mm = 1:length(pred)
          hlong = plot(50,pred(mm).perrmtx(1,ii),symb{mm}(1:2));
          set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
          hshort = plot(xval,pred(mm).perrmtx(2:end,ii),symb{mm});
          set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
        end

        if ii==1
          title('Local error (deg)','FontSize',kv.FontSize)
        end

        axis([5,75,26,54])
      end
      % Pooled
      axes(ha(2*ii+2));
      hlong = errorbar(50,mean(perrmtx(1,:)),std(perrmtx(1,:)),symbExp(1:2));
      set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)
      set(gca,'YAxisLocation','right','YMinorTick','on','TickLength',kv.TickLength)
      hold on
      hshort = errorbar(xval,mean(perrmtx(2:end,:),2),std(perrmtx(2:end,:),1,2),symbExp);
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      hlong = errorbar(50,mean(pred(mm).perrmtx(1,:),2),std(pred(mm).perrmtx(1,:),1,2),symb{mm}(1:2));
      set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
      hshort = errorbar(xval,mean(pred(mm).perrmtx(2:end,:),2),std(pred(mm).perrmtx(2:end,:),1,2),symb{mm});
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      axis([5,75,26,54])

      xlabel('SL (dB)','FontSize',kv.FontSize)


      % Quadrant error
      for ii=1:Nsub
        axes(ha(2*ii-1));

        hlong = plot(50,qerrmtx(1,ii),symbExp(1:2));

        set(gca,'YMinorTick','on')
    %     set(gca,'XTick',xval,'YTick',10:20:50)
        set(gca,'TickLength',kv.TickLength)

        set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)
        hold on
        hshort = plot(xval,qerrmtx(2:end,ii),symbExp);
        set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)

        for mm = 1:length(pred)
          hlong = plot(50,pred(mm).qerrmtx(1,ii),symb{mm}(1:2));
          set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
          hshort = plot(xval,pred(mm).qerrmtx(2:end,ii),symb{mm});
          set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
        end

        if ii==1
          title('Quadrant error (%)','FontSize',kv.FontSize)
        end

        axis([5,75,-4,54])

        % Listener ID
        ylabel([data(ii).id '            '],'FontSize',kv.FontSize,'Rotation',0,'FontWeight','bold');
    %     text(xtext,30,data(ii).id,'FontSize',kv.FontSize)
      end
      % Pooled
      axes(ha(2*ii+1));
      hlong = errorbar(50,mean(qerrmtx(1,:)),std(qerrmtx(1,:)),symbExp(1:2));
      set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)
      set(gca,'YMinorTick','on','TickLength',kv.TickLength)
      hold on
      hshort = errorbar(xval,mean(qerrmtx(2:end,:),2),std(qerrmtx(2:end,:),1,2),symbExp);
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      hlong = errorbar(50,mean(pred(mm).qerrmtx(1,:),2),std(pred(mm).qerrmtx(1,:),1,2),symb{mm}(1:2));
      set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
      hshort = errorbar(xval,mean(pred(mm).qerrmtx(2:end,:),2),std(pred(mm).qerrmtx(2:end,:),1,2),symb{mm});
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      axis([5,75,-4,54])
      ylabel('Pooled            ','FontSize',kv.FontSize,'Rotation',0,'FontWeight','bold');
      xlabel('SL (dB)','FontSize',kv.FontSize)

    end
    
    if flags.do_pmv
      
      conditions = data_baumgartner2015('ConditionNames');
      idcond = [1,2,5,8]; %{'Long','10dB','30dB','50dB','70dB'};
      Nc = length(idcond);
      
      hfig = figure;
      ha = tight_subplot(Nsub,Nc,kv.gap,kv.marg_h,kv.marg_w);
      
      for cc=1:Nc
        idc = idcond(cc);
        data = data_baumgartner2015(conditions{idc});
        
        for ll=1:Nsub
          axes( ha(cc + Nc*(ll-1)) )
          plot_baumgartner2014(pred.pmv{ll,idc}.p,pred.pmv{ll,idc}.tang,pred.pmv{ll,idc}.rang,...
              data(ll).itemlist(:,6),data(ll).itemlist(:,8),...
              'MarkerSize',kv.MarkerSize/2,'nocolorbar','cmax',0.05)
            
          % Labels
          set(gca,'FontSize',kv.FontSize-1)
          if not(cc==1)
            set(gca,'YTickLabel',[])
            ylabel('')
          elseif not(ll==4)
              ylabel('')
          end
          if not(ll==Nsub)
            set(gca,'XTickLabel',[])
            xlabel('')
          else
%             XTickLabel = get(gca,'XTickLabel');
%             XTickLabel(end:-2:1) = {' '};
            set(gca,'XTickLabelRotation',90)
            xlabel('')
          	if cc==3
              xlabel({' ';'Target Angle (deg)                                 '},'FontSize',kv.FontSize-1)
            end
          end
          if ll==1
            title(conditions{idc},'FontSize',kv.FontSize)
          end
          if cc==Nc
            set(gca,'YAxisLocation','right')
            ylabel(['            ' data(ll).id],'FontSize',kv.FontSize,'FontWeight','bold','Rotation',0);
          end
            
        end
        
      end
      
    end
  
  end
end

end



%% ------------------------------------------------------------------------
%  ---- INTERNAL FUNCTIONS ------------------------------------------------
%  ------------------------------------------------------------------------
function hM_warped = warp_hrtf(hM,fs)
% warps HRTFs acc. to Walder (2010)
% Usage: hM_warped = warp_hrtf(hM,fs)

N = fs;
fu = 2800;
fowarped = 8500;
fo = 16000;

fscala = [0:fs/N:fs-fs/N]';
hM_warped = zeros(512,size(hM,2),size(hM,3));
fuindex = max(find(fscala <= fu)); % 2800
fowindex = min(find(fscala >= fowarped));
foindex = min(find(fscala >= fo));

for canal = 1:size(hM,3)
   for el = 1:size(hM,2)
        yi = ones(fs/2+1,1)*(10^-(70/20));
        flin1 = [fscala(1:fuindex-1)];
        flin2 = [linspace(fscala(fuindex),fscala(foindex),fowindex-fuindex+1)]';
        fscalawarped = [flin1 ; flin2];

        % interpolate
        x = fscala(1:foindex);
        H = fft(hM(:,el,canal),N);
        Y = H(1:foindex);
        xi = fscalawarped;
        yi(1:length(xi),1) = interp1(x, Y, xi,'linear');

        yges=([yi; conj(flipud(yi(2:end-1)))]);
        hges = ifft([yges(1:end)],length(yges));
        hges=fftshift(ifft(yges));
        hwin=hges(fs/2-256:fs/2+768);
        hwinfade = FW_fade(hwin,512,24,96,192);
        hM_warped(1:end,el,canal)=hwinfade;
            
    end
end
end

function [syncrnfreq, GETtrain] = GETVocoder(filename,in,channum,lower,upper,alpha,GaussRate,stimpar)
warning('off')
% channel/timing parameters
srate=stimpar.SamplingRate;
nsamples=length(in); % length of sound record
duration=nsamples/srate; % duration of the signal (in s)
t=0:1/srate:duration;
t=t(1:nsamples);
extendedrange = 0;
if alpha == -1 % Log12ER case
    extendedrange = 1;
    alpha = 0.28;
elseif alpha == 0 % use predefined alphas
  switch channum
    case 3
      alpha=1.42;
    case 6
      alpha=0.67;
    case 9
      alpha=0.45;
    case 12
      alpha=0.33;
    case 18
      alpha=0.22;
    case 24
      alpha=0.17;
    otherwise 
      error(['Alpha not predefined for channels number of ' num2str(channu)]);
  end
end

% Synthesis: These are the crossover frequencies that the output signal is mapped to
crossoverfreqs = logspace( log10(lower), log10(upper), channum + 1);
if extendedrange == 1
    crossoverfreqs = [300,396,524,692,915,1209,1597,2110,2788,4200,6400,10000,16000];
end
syncrnfreq(:,1)=crossoverfreqs(1:end-1);
syncrnfreq(:,2)=crossoverfreqs(2:end);
for i=1:channum
    cf(i) = sqrt( syncrnfreq(i,1)*syncrnfreq(i,2) );
end

% Pulse train parameters
Gamma = alpha*cf;  
if extendedrange == 1
    Gamma(9) = 1412;
    Gamma(10) = 2200;
    Gamma(11) = 3600;
    Gamma(12) = 6000;
end
N = ceil(duration*GaussRate); % number of pulses
Genv=zeros(N,nsamples);
GETtrain=zeros(channum,nsamples);

% Make pulse trains
for i = 1:channum
    Teff = 1000/Gamma(i);
    if Teff > 3.75
    % if modulation depth is not 100%, make pulse train then modulate
        for n = 1:N
            % delay pulses by half a period so first Gaussian pulse doesn't
            % start at a max
            T = (n-0.5)/N*duration;
            Genv(n,:) = sqrt(Gamma(i)) * exp(-pi*(Gamma(i)*(t-T)).^2);
        end
        Genv_train(i,:) = sum(Genv);
        %modulate carrier
        GETtrain(i,:) = Genv_train(i,:) .* sin(2*pi*cf(i)*t);
        %normalize energy
        Energy(i) = norm(GETtrain(i,:))/sqrt(length(t));
%         Energy(i) = rms(GETtrain(i,:)); % !!!!!!!!!!!!
        GETtrain(i,:) = GETtrain(i,:)/Energy(i);
    else
    % if modulation depth is 100%, make modulated pulses and replicate
        T=(0.5)/N*duration;
        Genv=zeros(N,nsamples);
        Genv(1,:) = sqrt(Gamma(i)) * exp(-pi*(Gamma(i)*(t-T)).^2) .* sin(2*pi*cf(i)*t - T + pi/4); %!!! (t-T)
        Genv=repmat(Genv(1,:),[N 1]);
        for n=1:N
            T = round((n)/N*nsamples);
            Genv(n,:)=circshift(Genv(n,:),[1 T-1]);
        end
        GETtrain(i,:) = sum(Genv);
        %normalize energy
        Energy(i) = norm(GETtrain(i,:))/sqrt(length(t));
%         Energy(i) = rms(GETtrain(i,:)); % !!!!!!!!!!!!
        GETtrain(i,:) = GETtrain(i,:)/Energy(i);
    end
end

end

function out=channelize(fwavout, h, h0, in, channum, corners, syncrnfreq, ...
                        GETtrain, stimpar, amp, fadein, fadeout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** v1.2.0                                                           %
% GET Vocoder scaled by energy, not envelope                           %
%                                                                      %
% *** v1.1.0                                                           %
% Added Gaussian Envelope Tone (GET) Vocoder                           %
% GET pulse train is generated and passed to this function             %
% M. Goupell, May 2008                                                 %
%                                                                      %
% *** v1.0.0                                                           %
% Modified from ElecRang/matlab/makewav.m  v1.3.1                      %
% To be used with Loca by M. Goupell, Nov 2007                         % 
% Now program receives a sound rather than reading a speech file       %
% Rewrote according to the specifications (PM, Jan. 2008)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = hrtf
% h0 = hrtf for reference position
% in = reference noise
% noise = number of channels of noise vectors

srate=stimpar.SamplingRate;
N=length(in);               % length of sound record
d=0.5*srate;                % frequency scalar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read corner frequencies for channels                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis: These are the crossover frequencies that the output signal is mapped to
% calculated now in GETVocoder, passed to this function

% Analysis: These are the crossover frequencies that the input signal is subdivided by
if length(corners)<=channum
  error('You need at least one more corner frequency than number of channels');
end
anacrnfreq(:,1)=corners(1:end-1);
anacrnfreq(:,2)=corners(2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis: Filter Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inX=fftfilt(h,in);

order=4;     %order of butterworth filter
% out=zeros(1,N);
% filtX=zeros(channum,N);

for i=1:channum
   [b, a]=butter(order, anacrnfreq(i,:)/d);
   out=filter(b, a, inX);  % bandpass-filtering
%    E(i)=norm(out,1)/sqrt(N);
   E(i) = rms(out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% each channel
for i=1:channum
    outX(i,:) = GETtrain(i,:) * E(i);    
end
% sum me up scotty
out=sum(outX,1);

% for i = 1:channum
%     subplot(channum/3,3,i)
%     plot(outX(i,(1:480)))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save file                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out=out*(10^(amp/20))/sqrt(sum(out.^2))*sqrt(sum(in.^2))*sqrt(sum(h.^2))/sqrt(sum(h0.^2));
% disp(20*log10(sqrt(sum(out.^2))));
ii=max(max(abs(out)));
if ii>=1
  error(['Maximum amplitude value is ' num2str(20*log10(ii)) 'dB. Set the HRTF scaling factor lower to avoid clipping']);
end
out=FW_fade(out,0,fadein,fadeout);
% wavwrite(out,srate,stimpar.Resolution,fwavout);
end

function out = FW_fade(inp, len, fadein, fadeout, offset)
% FW_FADE crop/extend and fade in/out a vector.
%
% OUT = FW_FADE(INP, LEN, FADEIN, FADEOUT, OFFSET) crops or extends with zeros the signal INP
% up to length LEN. Additionally, the result is faded in/out using HANN window with 
% the length FADEIN/FADEOUT, respectively. If given, an offset can be added to show
% where the real signal begins.
%
% When used to crop signal, INP is cropped first, then faded out. 
% When used to extend signal, INP is faded out first, then extended too provide fading.
% 
% INP:     vector with signal (1xN or Nx1)
% LEN:     length of signal OUT (without OFFSET)
% FADEIN:  number of samples to fade in, beginning from OFFSET
% FADEOUT: number of samples to fade out, ending at the end of OUT (without OFFSET)
% OFFSET:  number of offset samples before FADEIN, optional
% OUT:     cropped/extended and faded vector
% 
% Setting a parameter to 0 disables corresponding functionality.

% ExpSuite - software framework for applications to perform experiments (related but not limited to psychoacoustics).
% Copyright (C) 2003-2010 Acoustics Research Institute - Austrian Academy of Sciences; Piotr Majdak and Michael Mihocic
% Licensed under the EUPL, Version 1.1 or ? as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence. 
% You may obtain a copy of the Licence at: http://ec.europa.eu/idabc/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% 7.11.2003
% 22.08.2005: improvement: INP may be 1xN or Nx1 now.
% Piotr Majdak (piotr@majdak.com)

ss=length(inp);    % get length of inp
	% offset
if ~exist('offset','var')
	offset = 0;
end
offset=round(offset);
len=round(len);
fadein=round(fadein);
fadeout=round(fadeout);
if offset>ss
	error('OFFSET is greater than signal length');
end
	% create new input signal discarding offset
inp2=inp(1+offset:end);
ss=length(inp2);

  % fade in
if fadein ~= 0 
  if fadein > ss
    error('FADEIN is greater than signal length');
  end
  han=hanning(fadein*2);
  if size(inp2,1)==1
    han=han';
  end
  inp2(1:fadein) = inp2(1:fadein).*han(1:fadein);
end
  % fade out window
  
if len == 0    
  len = ss;
end
if len <= ss
    % crop and fade out
  out=inp2(1:len);
    % fade out
  if fadeout ~= 0
    if fadeout > len
      error('FADEOUT is greater than cropped signal length');
    end
    han = hanning(2*fadeout);
    if size(out,1)==1
      han=han';
    end
    out(len-fadeout+1:end)=out(len-fadeout+1:end).*han(fadeout+1:end);
  end
else
    % fade out and extend
  if fadeout ~= 0
    if fadeout > ss
      error('FADEOUT is greater than signal length');
    end
    han = hanning(2*fadeout);
    inp2(ss-fadeout+1:end)=inp2(ss-fadeout+1:end).*han(fadeout+1:end);
  end
  out = [zeros(offset,1); inp2; zeros(len-ss,1)];
end
end

function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .1],[.1 .1])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
end

function s = gain2slope(g)
% s = gain2slope(g)

s = rad2deg(acos(1./sqrt(g.^2+1)));
end