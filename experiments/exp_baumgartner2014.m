function varargout = exp_baumgartner2014(varargin)
%EXP_BAUMGARTNER2014 Figures from Baumgartner et al. (2014)
%   Usage: data = exp_baumgartner2014(flag)
%
%   `exp_baumgartner2014(flag)` reproduces figures of the study from 
%   Baumgartner et al. (2014).
%
%   Optional fields of output *data* structure:
%
%   `data.contralateralGain`
%      contralateral gain of binaural weighting function
%
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'auto'     Re-calculate the file if it does not exist. Return 1 if the
%                file exist, otherwise 0. This is the default
%
%     'refresh'  Always recalculate the file.
%
%     'cached'   Always use the cached version. Throws an error if the
%                file does not exist.
%
%     'fig2'    Reproduce Fig.2:
%               Binaural weighting function best fitting results from 
%               Morimoto (2001) labeled as [1] and Macpherson and Sabin (2007) 
%               labeled as [2] in a squared error sense.
%
%     'fig3'    Reproduce Fig.3:
%               Partial and joint prediction residues, $e_\mathrm{PE}$ and/or 
%               $e_\mathrm{QE}$, as functions of the degree of selectivity, 
%               $\Gamma$, and the motoric response scatter, $\varepsilon$. 
%               Residuum functions are normalized to the minimum residuum 
%               obtained for the optimal parameter value. See text for details.
%               
%     'fig4'    Reproduce Fig.4:
%               Actual and predicted response patterns of listeners NH15, 
%               NH22, and NH72 (from left to right) when listening to  targets 
%               in the baseline condition and in the midsagittal plane. 
%               Actual response angles are shown as open circles. Probabilistic 
%               response predictions are encoded by brightness according to 
%               the color bar to the right. Actual (A) and predicted (P) 
%               performances (PEs and QEs) of the listener are listed above 
%               each panel.
%
%     'fig5'    Reproduce Fig.5:
%               Listener-specific baseline performance in the midsagittal plane. 
%               Local performance is shown in terms of PE (left), global 
%               performance in terms of QE (right). Correlation coefficients, 
%               $r$, and prediction residues, $e$, with respect to actual and 
%               predicted performances are shown above each panel.
%
%     'fig6'    Reproduce Fig.6:
%               Baseline performance as a function of the magnitudes of the 
%               lateral response angle. Symbols and whiskers show median values 
%               and inter-quartile ranges, respectively. Symbols were 
%               horizontally shifted to avoid overlaps. Correlation coefficients, 
%               $r$, and prediction residues, $e$, specify the correspondence 
%               between actual (A) and predicted (P) listener-specific performances. 
%               Predictions of the model without the sensory-motoric mapping (SMM) 
%               stage are shown with dashed lines.
%
%     'fig7'    Reproduce Fig.7:
%               Actual and predicted response patterns of listener NH62 when 
%               listening to broadband (left), low-pass filtered (center), 
%               or spectrally warped (right) DTFs of the midsagittal plane. 
%               Data were pooled within $\pm15^\circ$ of lateral angle.
%               All other conventions are as in Fig.4.
%
%     'fig8'    Reproduce Fig.8:
%               Effect of band limitation and spectral warping. 
%               Listeners were tested with broadband (BB), low-pass filtered (LP), 
%               and spectrally warped (W) DTFs. Actual: experimental results 
%               from Majdak et al. (2013). Part.: Model predictions for the 
%               actual eight participants based on the actually tested target 
%               positions. Pool: Model predictions for our pool of 18 listeners 
%               based on all possible target positions. Dotted horizontal lines 
%               represent chance rate. All other conventions are as in Fig.6.
%
%     'fig9'    Reproduce Fig.9:
%               Actual and predicted response patterns of exemplary listener 
%               NH12 listening to channel-limited HRTFs. Results for an unlimited 
%               number of channels (broadband click trains), 24, and 9 channels 
%               are shown from left to right. All other conventions are as in Fig.7.
%
%     'fig10'   Reproduce Fig.10:
%               Effect of spectral resolution in terms of varying the number 
%               of spectral channels used by a channel vocoder. Actual experimental 
%               results are from Goupell et al. (2010). Stimulation with broadband 
%               click trains (CL) represents an unlimited number of channels. 
%               All other conventions are as in Fig.8.  
%
%     'fig11'   Reproduce Fig.11:
%               Effect of non-individualized HRTFs, i.e., localizing with others' 
%               instead of own ears. Statistics summary replotted from Fig.13 
%               of Middlebrooks (1999b). Horizontal lines represent 25th, 50th, 
%               and 75th percentiles, the whiskers represent 5th and 95th percentiles, 
%               and crosses represent minima and maxima. Circles and squares 
%               represent mean values.
%
%     'fig12'   Reproduce Fig.12:
%               Effect of spectral ripples. Actual experimental results are 
%               from Macpherson et al. (2003). Either the ripple depth of 40dB (top) 
%               or the ripple density of one ripple/octave (bottom) was kept constant.
%               Predictions (P) of the model without the DCN stage are shown 
%               with dashed lines. All other conventions are as in Fig.8.  
%
%     'fig13'   Reproduce Fig.13:
%               Effect of high-frequency attenuation in speech. Actual experimental 
%               results are from Best et al. (2005). Open and filled symbols 
%               show actual and predicted results, respectively. Absolute polar 
%               angle errors (top) and QEs (bottom) averaged across listeners 
%               are shown.  
%
%   See also: baumgartner2013, data_baumgartner2013
%
%   Examples:
%   ---------
%
%   To display Fig.2 use :::
%
%     exp_baumgartner2014('fig2');
%
%   To display Fig.3 use :::
%
%     exp_baumgartner2014('fig3');
%
%   To display Fig.4 use :::
%
%     exp_baumgartner2014('fig4');
%
%   To display Fig.5 use :::
%
%     exp_baumgartner2014('fig5');
%
%   To display Fig.6 use :::
%
%     exp_baumgartner2014('fig6');
%
%   To display Fig.7 use :::
%
%     exp_baumgartner2014('fig7');
%
%   To display Fig.8 use :::
%
%     exp_baumgartner2014('fig8');
%
%   To display Fig.9 use :::
%
%     exp_baumgartner2014('fig9');
%
%   To display Fig.10 use :::
%
%     exp_baumgartner2014('fig10');
%
%   To display Fig.11 use :::
%
%     exp_baumgartner2014('fig11');
%
%   To display Fig.12 use :::
%
%     exp_baumgartner2014('fig12');
%
%   To display Fig.13 use :::
%
%     exp_baumgartner2014('fig13');
%
%   References: morimoto2001 macpherson2007

  
% AUTHOR: Robert Baumgartner


%% ------ Check input options --------------------------------------------

definput.import={'amtredofile'};
definput.flags.type = {'missingflag','fig2','fig3','fig4','fig5','fig6',...
                       'fig7','fig8','fig9','fig10','fig11','fig12','fig13'};
definput.flags.plot = {'plot','noplot'};


[flags]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

save_format='-v6';


%% General Plot Settings --------------------------------------------------

FontSize = 12;
MarkerSize = 6;

%% ------ FIG 2 -----------------------------------------------------------
if flags.do_fig2
  
  % Original Data
  Morimoto_left = [1,0.5,0];
  Morimoto_right = 1-Morimoto_left;
  Morimoto_ang = [60,0,-60];

  Macpherson_left = [3.75,0.5,-4]/10+0.5;
  Macpherson_right = [-3.75,-0.5,4.75]/10+0.5;
  Macpherson_ang = [30,0,-30];

  data =  [Morimoto_left, -Macpherson_right, Macpherson_left];
  lat =   [Morimoto_ang ,  Macpherson_ang  , Macpherson_ang];

  % Fit Slope
  bwslope = 1:0.1:90;
  resid = zeros(size(bwslope));
  for ii = 1:length(bwslope)
    binw_left = 1./(1+exp(-lat/bwslope(ii))); 
    resid(ii) = rms(binw_left-data);
  end
  [~,idmin] = min(resid);
  contralateralGain = bwslope(idmin);
  fprintf('Phi: %2.0f deg\n',contralateralGain)
  
  % Calculate specific weights to plot
  lat = -90:5:90;
  binw_left = 1./(1+exp(-lat/contralateralGain)); % weight of left ear signal with 0 <= binw <= 1
  binw_right = 1-binw_left;
  
  if flags.do_plot
    fig = figure;
    plot(lat,binw_left)
    hold on
    plot(lat,binw_right,'r--')
    plot(Morimoto_ang,Morimoto_left,'vk','MarkerSize',MarkerSize)
    plot(Macpherson_ang,Macpherson_left,'ok','MarkerSize',MarkerSize)

    plot(Morimoto_ang,Morimoto_left,'vb','MarkerSize',MarkerSize)
    plot(Morimoto_ang,Morimoto_right,'vr','MarkerSize',MarkerSize)
    plot(Macpherson_ang,Macpherson_left,'ob','MarkerSize',MarkerSize)
    plot(Macpherson_ang,Macpherson_right,'or','MarkerSize',MarkerSize)

    l = legend('\itL','\itR','[1]','[2]');
    set(l,'Location','East','FontSize',FontSize-1)
    set(gca,'XLim',[lat(1) lat(end)],'YLim',[-0.05 1.05],'XTick',-60:30:60,'FontSize',FontSize)
    xlabel('\phi_k (deg)','FontSize',FontSize)
    ylabel('w_{\zeta}(\phi_k)','FontSize',FontSize)
  end
  
  % Output
  data.contralateralGain = contralateralGain;
  
end

%% ------ FIG 3 -----------------------------------------------------------
if flags.do_fig3
  
  fn = [mfilename('fullpath'),'_parametrization.mat'];
  
  if amtredofile(fn,flags.redomode)
    
    autorefreshnotification(fn,flags)
    
    tempfn = fullfile(amtbasepath,'experiments','exp_baumgartner2014_parametrization'); % temporary folder
    mkdir(tempfn)
    
    gamma = [1,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,8,8,8,9,9,9,...
      10,10,10,12,12,12,16,16,16,30,30,30,100,100,100];
    mrs = [0,17,18,19,19,20,21,19,20,21, 0,10,100,17,19,20,21,23,30, 5,...
      19,20,21,19,20,21,19,20,21,19,20,21,19,20,21,19,20,21,19,20,21,19,20,21];

    for g = 1:length(gamma)

      latseg = -60:20:60; % centers of lateral segments
      dlat =  10;  % lateral range (+-) of each segment

      s = data_baumgartner2014('baseline','recalib','gamma',gamma(g),'mrsmsp',mrs(g));

      qe_exp = zeros(length(s),length(latseg));
      pe_exp = zeros(length(s),length(latseg));
      for ll = 1:length(s)

        s(ll).target = [];
        s(ll).response = [];
        s(ll).Nt = [];
        for ii = 1:length(latseg)

          latresp = s(ll).mm1(:,7);
          idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
          s(ll).mm2 = s(ll).mm1(idlat,:);

          s(ll).mm2(:,7) = 0; % set lateral angle to 0deg such that localizationerror works also outside +-30deg

          pe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'rmsPmedianlocal'));
          qe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'querrMiddlebrooks'));

          s(ll).target{ii} = real(s(ll).mm2(:,6)); % polar angle of target
          s(ll).response{ii} = real(s(ll).mm2(:,8)); % polar angle of response
          s(ll).Nt{ii} = length(s(ll).target{ii});

        end
      end


      %% LocaMo
      qe = zeros(length(s),length(latseg));
      pe = zeros(length(s),length(latseg));
      for ll = 1:length(s)

        for ii = 1:length(latseg)

          s(ll).sphrtfs{ii} = 0;     % init
          s(ll).p{ii} = 0;        % init

          [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
          [s(ll).p{ii},respangs] = baumgartner2014(...
              s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
              'S',s(ll).S,'lat',latseg(ii),'polsamp',polang,...
              'gamma',gamma(g),'mrsmsp',mrs(g)); 

          if s(ll).Nt{ii} > 0
            [ qe(ll,ii),pe(ll,ii) ] = pmv2ppp( ...
                s(ll).p{ii} , polang , respangs , s(ll).target{ii});
          else
            qe(ll,ii) = NaN; 
            pe(ll,ii) = NaN;
          end

        end

      end
      s = rmfield(s,{'Obj','mm1','mm2','sphrtfs'}); % reduce file size
      savename = ['result_baseline_g' num2str(gamma(g),'%u') '_mrs' num2str(mrs(g),'%u')];
      save(fullfile(tempfn,savename),'s', 'qe', 'pe', 'qe_exp', 'pe_exp', 'latseg')
    end

    %% Combine results to single mat file
    fname = dir(fullfile(tempfn,'result_baseline_g*.mat'));
    perr = zeros(18,7,length(fn));
    qerr = perr;
    for ii = 1:length(fname)
      tmp = load(fullfile(tempfn,fname(ii).name));
      perr(:,:,ii) = tmp.pe;
      qerr(:,:,ii) = tmp.qe;

      if ii == 1
        perr_exp = tmp.pe_exp;
        qerr_exp = tmp.qe_exp;
      end

      gKey = '_g';
      gIndex = strfind(fname(ii).name,gKey);
      g(ii) = sscanf(fname(ii).name(gIndex(1) + length(gKey):end), '%g', 1);

      mrsKey = '_mrs';
      mrsIndex = strfind(fname(ii).name,mrsKey);
      mrs(ii) = sscanf(fname(ii).name(mrsIndex(1) + length(mrsKey):end), '%g', 1);
    end
    
    % Sort data acc. to ascending gamma
    [g,ig] = sort(g);
    mrs = mrs(ig);
    perr = perr(:,:,ig);
    qerr = qerr(:,:,ig);

    % Number of targets for each listeners (1st dim) and lateral segment (2nd
    % dim)

    Ntargets = zeros(18,7);
    for jj = 1:18
      Ntargets(jj,:) = [tmp.s(jj).Nt{:}];
    end
    
    save(fn,'perr','perr_exp','qerr','qerr_exp','g','mrs','Ntargets',save_format);
    
    rmdir(tempfn,'s')
    
  else
    load(fn);
  end
  varargout{1} = load(fn);
  
  [qerr0,perr0] = pmv2ppp(ones(49,44)); % chance performances

  % extract all different gammas
  gamma = g(1);
  for ii =2:length(g)
    if g(ii) > gamma(end)
      gamma = [gamma;g(ii)];
    end
  end

  Nlat = size(perr,2);
  Nset = size(perr,3);
  idnum = Ntargets ~= 0;
  relNt = Ntargets/sum(Ntargets(:));


  % Compute all residues
  resid.perr = zeros(length(g),1);
  resid.qerr = resid.perr;
  for ii = 1:Nset

    dperr = perr_exp - perr(:,:,ii);
    dqerr = qerr_exp - qerr(:,:,ii);
    resid.perr(ii) = sqrt( relNt(idnum)' * (dperr(idnum)).^2 );
    resid.qerr(ii) = sqrt( relNt(idnum)' * (dqerr(idnum)).^2 );

  end
  resid.total = resid.perr/perr0 + resid.qerr/qerr0;

  % Select optimal residues for various gamma
  id_g = zeros(length(gamma),1);
  etotal_g = zeros(length(gamma),1);
  for ii = 1:length(gamma)
    idgamma = find(g == gamma(ii));
    [etotal_g(ii),id] = min(resid.total(idgamma));
    id_g(ii) = idgamma(id);
  end
  eperr_g = resid.perr(id_g);
  eqerr_g = resid.qerr(id_g);

  [tmp,idopt] = min(etotal_g);
  etotal_g = etotal_g / etotal_g(idopt);
  eperr_g = eperr_g / eperr_g(idopt);
  eqerr_g = eqerr_g / eqerr_g(idopt);


  % Select residues for optimal gamma and various mrs
  idgammaopt = find(g == gamma(idopt));
  mrs_gopt = mrs(idgammaopt);
  [mrssort,idsort] = sort(mrs_gopt);
  idmrs = idgammaopt(idsort);
  [tmp,idopt_mrs] = min(resid.total(idmrs));
  idnorm = idmrs(idopt_mrs);
  etotal_gopt = resid.total(idmrs) / resid.total(idnorm);
  eperr_gopt = resid.perr(idmrs) / resid.perr(idnorm);
  eqerr_gopt = resid.qerr(idmrs) / resid.qerr(idnorm);
  
  
  if flags.do_plot
    
    %% Plot residues for various gamma

    % Interpolate data
    gamma_int = logspace(0,2.1,1000);
    inttype = 'cubic';
    dperr_int = interp1(log10(gamma),eperr_g,log10(gamma_int),inttype);
    dqerr_int = interp1(log10(gamma),eqerr_g,log10(gamma_int),inttype);
    dtot_int = interp1(log10(gamma),etotal_g,log10(gamma_int),inttype);

    % Plot
    fig=figure;
    subplot(1,2,1)
    semilogx(gamma_int,dperr_int,'k: ')
    hold on
    semilogx(gamma_int,dqerr_int,'k--')
    semilogx(gamma_int,dtot_int,'k-')
    semilogx(gamma(idopt),0.95,'vk','MarkerFaceColor','k','MarkerSize',MarkerSize+1)

    leg = legend('PE','QE','PE+QE','\{\epsilon,\Gamma\}_{opt}');
    set(leg,'Location','northeast','FontSize',FontSize-1)

    ylabel('e(\Gamma) / e(\Gamma_{opt})','FontSize',FontSize)
    xlabel('\Gamma','FontSize',FontSize)

    set(gca,'XLim',[gamma(1)-0.1 gamma(end)+20],'YLim',[0.91 1.7],'XMinorTick','on',...
      'FontSize',FontSize-1,'XTickLabel',[1 10 100])

    %% Plot residues for optimal gamma and various mrs

    % Interpolation
    mrs_int = 0:0.1:45;
    inttype = 'cubic';
    dperr_int = interp1(mrssort,eperr_gopt,mrs_int,inttype);
    dqerr_int = interp1(mrssort,eqerr_gopt,mrs_int,inttype);
    dtot_int = interp1(mrssort,etotal_gopt,mrs_int,inttype);

    % Plot
    subplot(1,2,2)
    plot(mrs_int,dperr_int,'k:')
    hold on
    plot(mrs_int,dqerr_int,'k--')
    plot(mrs_int,dtot_int,'k-')
    plot(mrs(idopt_mrs),0.95,'vk','MarkerFaceColor','k','MarkerSize',MarkerSize+1)

    ylabel('e(\epsilon) / e(\epsilon_{opt})','FontSize',FontSize)
    xlabel('\epsilon (\circ)','FontSize',FontSize)

    set(gca,'XLim',[mrssort(1) 40],'YLim',[0.91 1.7],'XMinorTick','on',...
      'FontSize',FontSize-1)
    
  end
end

%% ------ FIG 4 -----------------------------------------------------------
if flags.do_fig4
  
  latseg = 0; % centers of lateral segments
  dlat =  10;  % lateral range (+-) of each segment

  s = data_baumgartner2014('baseline');
  
  idselect = ismember({s.id},{'NH15','NH22','NH72'});
  s = s(idselect);

  qe_exp = zeros(length(s),length(latseg));
  pe_exp = zeros(length(s),length(latseg));
  
  for ll = 1:length(s)

    s(ll).target = [];
    s(ll).response = [];
    s(ll).Nt = [];
    for ii = 1:length(latseg)
      
      latresp = s(ll).mm1(:,7);
      idlat = latresp <= latseg(ii)+dlat & latresp >= latseg(ii)-dlat;
      s(ll).mm2 = s(ll).mm1(idlat,:);
      
      % set lateral angle to 0deg such that localizationerror works outside +-30deg
      s(ll).mm2(:,7) = 0; 

      pe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'rmsPmedianlocal'));
      qe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'querrMiddlebrooks'));

      s(ll).target{ii} = real(s(ll).mm2(:,6)); % polar angle of target
      s(ll).response{ii} = real(s(ll).mm2(:,8)); % polar angle of response
      s(ll).Nt{ii} = length(s(ll).target{ii});

    end
  end


  %% LocaMo
  qe = zeros(length(s),length(latseg));
  pe = zeros(length(s),length(latseg));
  for ll = 1:length(s)

    for ii = 1:length(latseg)

      s(ll).sphrtfs{ii} = 0;     % init
      s(ll).p{ii} = 0;        % init

      [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
      [s(ll).p{ii},respangs] = baumgartner2014(...
          s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
          'S',s(ll).S,'lat',latseg(ii),'polsamp',polang); 

      if s(ll).Nt{ii} > 0
        [ qe(ll,ii),pe(ll,ii) ] = pmv2ppp( ...
            s(ll).p{ii} , polang , respangs , s(ll).target{ii});
      else
        qe(ll,ii) = NaN; 
        pe(ll,ii) = NaN;
      end

      if flags.do_plot
        if ll ==1; figure; end
        subplot(1,3,ll)
        idend = min(150,length(s(ll).target{ii}));
        plotbaumgartner2013(s(ll).p{ii},polang,respangs,...
                  s(ll).target{ii}(1:idend),s(ll).response{ii}(1:idend),...
                  'MarkerSize',MarkerSize,'cmax',0.05,'nocolorbar');
        title({['A: PE = ' num2str(pe_exp(ll,ii),2) 'deg, QE = ' num2str(qe_exp(ll,ii),2) '%'];['P: PE = ' num2str(pe(ll,ii),2) 'deg, QE = ' num2str(qe(ll,ii),2) '%']},...
          'FontSize',FontSize-1)
        xlabel('Target Angle (deg)','FontSize',FontSize)
        ylabel('Response Angle (deg)','FontSize',FontSize)
        set(gca,'FontSize',FontSize-1)
        set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
        set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})
      end
    end

  end
  
  s = rmfield(s,{'Obj','mm1','mm2','sphrtfs'}); % reduce file size 
  
  varargout{1} = s;
  
end

%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5
  
  fn = [mfilename('fullpath'),'_baseline.mat'];
  
  if amtredofile(fn,flags.redomode)
    autorefreshnotification(fn,flags)
    
    latseg = -60:20:60; % centers of lateral segments
    dlat =  10;  % lateral range (+-) of each segment

    s = data_baumgartner2014('baseline');

    qe_exp = zeros(length(s),length(latseg));
    pe_exp = zeros(length(s),length(latseg));
    for ll = 1:length(s)

      s(ll).target = [];
      s(ll).response = [];
      s(ll).Nt = [];
      for ii = 1:length(latseg)
        
        latresp = s(ll).mm1(:,7);
        idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
        s(ll).mm2 = s(ll).mm1(idlat,:);

        s(ll).mm2(:,7) = 0; % set lateral angle to 0deg such that localizationerror works outside +-30deg

        pe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'rmsPmedianlocal'));
        qe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'querrMiddlebrooks'));

        s(ll).target{ii} = real(s(ll).mm2(:,6)); % polar angle of target
        s(ll).response{ii} = real(s(ll).mm2(:,8)); % polar angle of response
        s(ll).Nt{ii} = length(s(ll).target{ii});

      end
    end


    %% LocaMo
    qe = zeros(length(s),length(latseg));
    pe = zeros(length(s),length(latseg));
    for ll = 1:length(s)

      for ii = 1:length(latseg)

        s(ll).sphrtfs{ii} = 0;     % init
        s(ll).p{ii} = 0;        % init

        [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
        [s(ll).p{ii},respangs] = baumgartner2014(...
            s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
            'S',s(ll).S,'lat',latseg(ii),'polsamp',polang); 

        if s(ll).Nt{ii} > 0
          [ qe(ll,ii),pe(ll,ii) ] = pmv2ppp( ...
              s(ll).p{ii} , polang , respangs , s(ll).target{ii});
        else
          qe(ll,ii) = NaN; 
          pe(ll,ii) = NaN;
        end

      end

    end

    s = rmfield(s,{'Obj','mm1','mm2','sphrtfs'}); % reduce file size 
    save(fn,'s', 'qe', 'pe', 'qe_exp', 'pe_exp', 'latseg',save_format);
    
  end
  
  model{1} = load(fn);
  varargout{1} = load(fn);
  
  sign = {'ko';'kd';'kd';'k<';'kp';'k>'};

  latseg = model{1}.latseg; 
  idlat = latseg == 0;
  pe_exp = model{1}.pe_exp(:,idlat);
  qe_exp = model{1}.qe_exp(:,idlat);
  [tmp,idsort] = sort(pe_exp);

  Ns = length(model{1}.s);
  for ll = 1:Ns
      NHs(ll,:) = model{1}.s(ll).id;
  end

  % Mean RMS Differences
  dpe = zeros(1,length(model));
  dqe = dpe;
  for ii = 1:length(model)
    idnum = not(isnan(model{ii}.pe_exp(idsort,idlat)) | isnan(model{ii}.pe(idsort,idlat)));
    idsort = idsort(idnum);
    Ns = length(model{ii}.s);
    relfreq = zeros(Ns,1);
    Ntall = zeros(Ns,1);
    for jj = 1:Ns
      Ntlat = [model{ii}.s(jj).Nt{idlat}];
      Ntall(jj) = sum(Ntlat);
      relfreq(jj,1) = Ntlat/Ntall(jj);
    end
    relfreq = relfreq.*repmat(Ntall,1,1)/sum(Ntall);
    dpe(ii) = sqrt( relfreq(idsort)' * (model{ii}.pe_exp(idsort,idlat) - model{ii}.pe(idsort,idlat)).^2 );
    dqe(ii) = sqrt( relfreq(idsort)' * (model{ii}.qe_exp(idsort,idlat) - model{ii}.qe(idsort,idlat)).^2 );
  end

  % Correlation Coefficients
  for ii = 1:length(model)
    idnum = not(isnan(model{ii}.pe_exp(:,idlat)) | isnan(model{ii}.pe(:,idlat)));
    r = corrcoef(model{ii}.pe_exp(idnum,idlat),model{ii}.pe(idnum,idlat));
    r_pe(ii) = r(2);
    idnum = not(isnan(model{ii}.qe_exp(:,idlat)) | isnan(model{ii}.qe(:,idlat)));
    [r,p] = corrcoef(model{ii}.qe_exp(idnum,idlat),model{ii}.qe(idnum,idlat));
    r_qe(ii) = r(2);
  end

  if flags.do_plot
    fig=figure;
    
    %% Polar Error

    subplot(1,2,1)

    ylim = [23.1 43];

    h = bar(pe_exp(idsort));
    hold on
    for ii = 1:length(model)
      plot(model{ii}.pe(idsort,idlat),sign{ii},'MarkerSize',MarkerSize,'MarkerFaceColor','k')
    end
    set(gca,'XLim',[0 Ns+1],'XTick',1:Ns,'XTickLabel',NHs(idsort,3:4),...
        'YLim',ylim,'YMinorTick','on','FontSize',FontSize-2)
      plot([0,Ns+1],[ylim(1),ylim(1)]+0.03,'k')
    set(h,'FaceColor','white')
    xlabel('Listener (NH)','FontSize',FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',FontSize)

      l = legend('Actual','Predicted');
      set(l,'FontSize',FontSize-2,'Location','northwest')
      title(['e_{PE} = ' num2str(dpe,'%0.1f') '\circ ; r_{PE} = ' num2str(r_pe,'%0.2f')],...
      'FontSize',FontSize)

    %% Quadrant Error

    subplot(1,2,2)

    ylim = [0 24.5];

    h = bar(qe_exp(idsort));
    hold on
    for ii = 1:length(model)
      plot(model{ii}.qe(idsort,idlat),sign{ii},'MarkerSize',MarkerSize,'MarkerFaceColor','k')
    end

    set(gca,'XLim',[0 Ns+1],'XTick',1:Ns,'XTickLabel',NHs(idsort,3:4),...
        'YLim',ylim,'YMinorTick','on','FontSize',FontSize-2)
    set(h,'FaceColor','white')
    xlabel('Listener (NH)','FontSize',FontSize)
    ylabel('Quadrant Error (%)','FontSize',FontSize)

    title(['e_{QE} = ' num2str(dqe,'%0.1f') '% ; r_{QE} = ' num2str(r_qe,'%0.2f')],...
      'FontSize',FontSize)

    set(fig,'PaperPosition',[1,1,9,3.5])
  end

end

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6
  
  fn = [mfilename('fullpath'),'_baseline.mat'];
  
  if amtredofile(fn,flags.redomode)
    autorefreshnotification(fn,flags)
    if flags.do_refresh
      exp_baumgartner2014('fig5','refresh');
    else
      exp_baumgartner2014('fig5');
    end
  else
    varargout{1} = load(fn);
    load(fn);
  end
  
  paradata = load([mfilename('fullpath'),'_parametrization.mat']);
  idmrs0 = paradata.g == 6 & paradata.mrs == 0;
  mrs0.pe = paradata.perr(:,:,idmrs0);
  mrs0.qe = paradata.qerr(:,:,idmrs0);

  dx = 2;

  %% # of targets
  Ns = length(s);
  Nlat = length(latseg);
  Ntlat = zeros(Ns,Nlat);
  relfreq = zeros(Ns,Nlat);
  Ntall = zeros(Ns,1);
  for jj = 1:Ns
    Ntlat(jj,:) = [s(jj).Nt{:}];
    Ntall(jj) = sum(Ntlat(jj,:));
    relfreq(jj,:) = Ntlat(jj,:)/Ntall(jj);
  end
  relfreq = relfreq.*repmat(Ntall,1,Nlat)/sum(Ntall);

  %% Pooling to lateralization
  idlat0 = round(Nlat/2);
  idleft = idlat0-1:-1:1;
  idright = idlat0+1:Nlat;
  latseg = latseg(idlat0:end);
  relfreqLR = Ntlat(:,idleft) ./ (Ntlat(:,idleft) + Ntlat(:,idright) + eps);

  pe = [pe(:,idlat0) , relfreqLR.*pe(:,idleft) + (1-relfreqLR).*pe(:,idright)];
  pe_exp = [pe_exp(:,idlat0) , relfreqLR.*pe_exp(:,idleft) + (1-relfreqLR).*pe_exp(:,idright)];
  qe = [qe(:,idlat0) , relfreqLR.*qe(:,idleft) + (1-relfreqLR).*qe(:,idright)];
  qe_exp = [qe_exp(:,idlat0) , relfreqLR.*qe_exp(:,idleft) + (1-relfreqLR).*qe_exp(:,idright)];
  relfreq = [relfreq(:,idlat0) , relfreq(:,1:idlat0-1) + relfreq(:,Nlat:-1:idlat0+1)];

  mrs0.pe = [mrs0.pe(:,idlat0) , relfreqLR.*mrs0.pe(:,idleft) + (1-relfreqLR).*mrs0.pe(:,idright)];
  mrs0.qe = [mrs0.qe(:,idlat0) , relfreqLR.*mrs0.qe(:,idleft) + (1-relfreqLR).*mrs0.qe(:,idright)];

  %% Evaluation Metrics
  idnum = not(isnan(pe_exp) | isnan(pe));
  dpe = sqrt( relfreq(idnum)' * (pe_exp(idnum) - pe(idnum)).^2 );
  dqe = sqrt( relfreq(idnum)' * (qe_exp(idnum) - qe(idnum)).^2 );
  [r_pe,p_pe] = corrcoef(pe_exp(idnum),pe(idnum));
  [r_qe,p_qe] = corrcoef(qe_exp(idnum),qe(idnum));

  mrs0.dpe = sqrt( relfreq(idnum)' * (pe_exp(idnum) - mrs0.pe(idnum)).^2 );
  mrs0.dqe = sqrt( relfreq(idnum)' * (qe_exp(idnum) - mrs0.qe(idnum)).^2 );
  [mrs0.r_pe,mrs0.p_pe] = corrcoef(pe_exp(idnum),mrs0.pe(idnum));
  [mrs0.r_qe,mrs0.p_qe] = corrcoef(qe_exp(idnum),mrs0.qe(idnum));


  %% Quartiles
  quart_pe = zeros(3,length(latseg),2); % 1st dim: 25/50/75 quantiles; 2nd dim: lat; 3rd dim: model/experiment/mrs0
  quart_qe = zeros(3,length(latseg),2);
  qlow = 0.25;
  qhigh = 0.75;
  for ii = 1:length(latseg)

    id = not(isnan(pe(:,ii)));
    quart_pe(:,ii,1) = quantile(pe(id,ii),[qlow .50 qhigh]);
    quart_pe(:,ii,3) = quantile(mrs0.pe(id,ii),[qlow .50 qhigh]);
    id = not(isnan(pe_exp(:,ii)));
    quart_pe(:,ii,2) = quantile(pe_exp(id,ii),[qlow .50 qhigh]);

    id = not(isnan(qe(:,ii)));
    quart_qe(:,ii,1) = quantile(qe(id,ii),[qlow .50 qhigh]);
    quart_qe(:,ii,3) = quantile(mrs0.qe(id,ii),[qlow .50 qhigh]);
    id = not(isnan(qe_exp(:,ii)));
    quart_qe(:,ii,2) = quantile(qe_exp(id,ii),[qlow .50 qhigh]);

  end


  if flags.do_plot
    
    %% PE

    fig = figure;
    subplot(1,2,1)
    errorbar(latseg-dx,quart_pe(2,:,1),...
      quart_pe(2,:,1)-quart_pe(1,:,1),...
      quart_pe(3,:,1)-quart_pe(2,:,1),...
      'ok-','MarkerSize',MarkerSize,'MarkerFaceColor','k');
    hold on
    errorbar(latseg,quart_pe(2,:,3),...
      quart_pe(2,:,3)-quart_pe(1,:,3),...
      quart_pe(3,:,3)-quart_pe(2,:,3),...
      'dk--','MarkerSize',MarkerSize-1,'MarkerFaceColor','k');
    errorbar(latseg+dx,quart_pe(2,:,2),...
      quart_pe(2,:,2)-quart_pe(1,:,2),...
      quart_pe(3,:,2)-quart_pe(2,:,2),...
      'ok-','MarkerSize',MarkerSize,'MarkerFaceColor','w');

    title(['e_{PE} = ' num2str(dpe,'%0.1f') '\circ ; r_{PE} = ' num2str(r_pe(2),'%0.2f')],...
      'FontSize',FontSize)
    set(gca,'XLim',[min(latseg)-2*dx,max(latseg)+2*dx],'YLim',[22.1,45.9],...
      'YMinorTick','on','FontSize',FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',FontSize)
    xlabel('Magnitude of Lateral Angle (deg)','FontSize',FontSize)


    %% QE

    subplot(1,2,2)
    errorbar(latseg-dx,quart_qe(2,:,1),...
      quart_qe(2,:,1)-quart_qe(1,:,1),...
      quart_qe(3,:,1)-quart_qe(2,:,1),...
      'ok-','MarkerSize',MarkerSize,'MarkerFaceColor','k');
    hold on
    errorbar(latseg,quart_qe(2,:,3),...
      quart_qe(2,:,3)-quart_qe(1,:,3),...
      quart_qe(3,:,3)-quart_qe(2,:,3),...
      'dk--','MarkerSize',MarkerSize-1,'MarkerFaceColor','k');
    errorbar(latseg+dx,quart_qe(2,:,2),...
      quart_qe(2,:,2)-quart_qe(1,:,2),...
      quart_qe(3,:,2)-quart_qe(2,:,2),...
      'ok-','MarkerSize',MarkerSize,'MarkerFaceColor','w');
    title(['e_{QE} = ' num2str(dqe,'%0.1f') '% ; r_{QE} = ' num2str(r_qe(2),'%0.2f')],...
      'FontSize',FontSize)
    set(gca,'XLim',[min(latseg)-2*dx,max(latseg)+2*dx],'YLim',[2.1,27.5],...
      'XTick',latseg,'YMinorTick','on','FontSize',FontSize)
    ylabel('Quadrant Error (%)','FontSize',FontSize)
    xlabel('Magnitude of Lateral Angle (deg)','FontSize',FontSize)

    l = legend('P (with SMM)','P (w/o SMM)','A');
    set(l,'FontSize',FontSize-1,'Location','northwest')

    set(fig,'PaperPosition',[1,1,10,3.5])
    
  end
end

%% ------ FIG 7 -----------------------------------------------------------
if flags.do_fig7
  
  latdivision = 0;  % lateral angle
  dlat = 15;

  % Experimental Settings
  Conditions = {'BB','LP','W'};


  %% Computations
  s = data_baumgartner2014('pool');  
  s = s(14); % for PMV plot of NH62
  chance = [];
  for C = 1:length(Conditions)

    Cond = Conditions{C};

    %% Data

    % Experimental data
    data = data_majdak2013(Cond);
    for ll = 1:length(s)
      if sum(ismember({data.id},s(ll).id)) % participant ?
        s(ll).mm1=data(ismember({data.id},s(ll).id)).mtx;
        for ii = 1:length(latdivision)
          latresp = s(ll).mm1(:,7);
          idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
          mm2 = s(ll).mm1(idlat,:);
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
                load exp_baumgartner2014_spatstrat_lpfilter
                s(ll).spdtfs_c{ii} = filter(blp,alp,s(ll).spdtfs{ii});
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

        [s(ll).p{ii},rang] = baumgartner2014(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
              'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii});

        [ qe(ii),pe(ii) ] = pmv2ppp(s(ll).p{ii} , s(ll).polang{ii} , rang);

        if sum(ismember({data.id},s(ll).id)) % if participant then actual targets
          [ qe_t(ii),pe_t(ii) ] = pmv2ppp( ...
              s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

        end

      end

      % Model results of pool
      s(ll).qe_pool(C,1) = mean(qe); 
      s(ll).pe_pool(C,1) = mean(pe);

      if sum(ismember({data.id},s(ll).id)) % participant ?
        % Actual experimental results
        s(ll).qe_exp(C,1) = localizationerror(s(ll).mm1,'querrMiddlebrooks');
        s(ll).pe_exp(C,1) = localizationerror(s(ll).mm1,'rmsPmedianlocal');
        s(ll).Nt(C,1) = size(s(ll).mm1,1);
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
          plotbaumgartner2013(s(ll).p{ii},s(ll).polang{ii},rang,...
                targets,responses,'MarkerSize',MarkerSize,'cmax',0.05,'nocolorbar')
          Nt = length(targets);
          tmp.m = [zeros(Nt,5) targets(:) zeros(Nt,1) responses(:)];
          tmp.qe = localizationerror(tmp.m,'querrMiddlebrooks');
          tmp.pe = localizationerror(tmp.m,'rmsPmedianlocal');
          title({['A: PE = ' num2str(tmp.pe,2) 'deg, QE = ' num2str(tmp.qe,2) '%'];...
            ['P: PE = ' num2str(pe_t(ii),2) 'deg, QE = ' num2str(qe_t(ii),2) '%']},'FontSize',FontSize-1)
          xlabel('Target Angle (deg)','FontSize',FontSize)
          ylabel('Response Angle (deg)','FontSize',FontSize)
          set(gca,'FontSize',FontSize-1)
          set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
          set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})

        end

      end

    end

  end
  
  varargout{1} = s;

end


%% ------ FIG 8 -----------------------------------------------------------
if flags.do_fig8
  
  fn = [mfilename('fullpath'),'_spatstrat.mat'];
  
  if amtredofile(fn,flags.redomode)
    
    autorefreshnotification(fn,flags)
    
    latdivision = [-20,0,20];            % lateral angle
    dlat = 10;

    % Experimental Settings
    Conditions = {'BB','LP','W'};

    %% Computations
    s = data_baumgartner2014('pool');
    chance = [];
    for C = 1:length(Conditions)

      Cond = Conditions{C};

      %% Data

      % Experimental data
      data = data_majdak2013(Cond);
      for ll = 1:length(s)
        if sum(ismember({data.id},s(ll).id)) % if actual participant
          s(ll).mm1=data(ismember({data.id},s(ll).id)).mtx; 
          for ii = 1:length(latdivision)
            latresp = s(ll).mm1(:,7);
            idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
            mm2 = s(ll).mm1(idlat,:);
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
                fn_filters = 'exp_baumgartner2014_spatstrat_lpfilter.mat';
                if not(exist(fn_filters,'file'))
                  disp(['Downloading ' fn_filters ' from http://www.kfs.oeaw.ac.at/']);
                  targetfn = fullfile(amtbasepath,'humandata',fn_filters);
                  sourcefn = ['http://www.kfs.oeaw.ac.at/research/experimental_audiology/projects/amt/' fn_filters];
                  urlwrite(sourcefn,targetfn);
                end
                load(fn_filters)
                s(ll).spdtfs_c{ii} = filter(blp,alp,s(ll).spdtfs{ii});
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

          [s(ll).p{ii},rang] = baumgartner2014(...
                s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
                'S',s(ll).S,'lat',latdivision(ii),...
                'polsamp',s(ll).polang{ii});

          [ qe(ii),pe(ii) ] = pmv2ppp(s(ll).p{ii} , s(ll).polang{ii} , rang);

          if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
            [ qe_t(ii),pe_t(ii) ] = pmv2ppp( ...
                s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

          end

        end

        % Model results of pool
        wlat = cos(deg2rad(latdivision)); % weighting compensating lateral compression
        wlat = wlat/sum(wlat);
        s(ll).qe_pool(C,1) = wlat * qe(:); 
        s(ll).pe_pool(C,1) = wlat * pe(:);

        if sum(ismember({data.id},s(ll).id)) % if actual participant 
          % Actual experimental results
          s(ll).qe_exp(C,1) = localizationerror(s(ll).mm1,'querrMiddlebrooks');
          s(ll).pe_exp(C,1) = localizationerror(s(ll).mm1,'rmsPmedianlocal');
          s(ll).Nt(C,1) = size(s(ll).mm1,1);
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

    %% Compute Chance Performance
    chance = repmat(chance,10,1);
    id_chance = randi(size(chance,1),size(chance,1),1);
    chance(:,8) = chance(id_chance,6);
    pe_chance = localizationerror(chance,'rmsPmedianlocal');
    qe_chance = localizationerror(chance,'querrMiddlebrooks');

    [r,p] =  corrcoef([s.qe_exp],[s.qe_part]);
    cc.qe.r = r(2);
    cc.qe.p = p(2);
    fprintf('QE: r = %0.2f, p = %0.3f\n',r(2),p(2));

    [r,p] =  corrcoef([s.pe_exp],[s.pe_part]);
    cc.pe.r = r(2);
    cc.pe.p = p(2);
    fprintf('PE: r = %0.2f, p = %0.3f\n',r(2),p(2));

    save(fn,'cc', 's', 'pe_chance', 'qe_chance',save_format);
  else
    load(fn);
  end
  varargout{1} = load(fn);
  
  %% Measures

  % Quartiles
  quart_pe_part = quantile([s.pe_part]',[.25 .50 .75]);
  quart_qe_part = quantile([s.qe_part]',[.25 .50 .75]);

  quart_pe_pool = quantile([s.pe_pool]',[.25 .50 .75]);
  quart_qe_pool = quantile([s.qe_pool]',[.25 .50 .75]);

  quart_pe_exp = quantile([s.pe_exp]',[.25 .50 .75]);
  quart_qe_exp = quantile([s.qe_exp]',[.25 .50 .75]);

  % RMS Differences
  % individual:
  Ntargets = [s.Nt]'; % # of targets
  relfreq = Ntargets/sum(Ntargets(:));
  sd_pe = ([s.pe_part]'-[s.pe_exp]').^2; % squared differences
  dpe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
  sd_qe = ([s.qe_part]'-[s.qe_exp]').^2;
  dqe = sqrt(relfreq(:)' * sd_qe(:));

  % Chance performance
  qe0 = qe_chance;
  pe0 = pe_chance;

  if flags.do_plot
    
    dx = 0.1;
    
    figure 

    subplot(121)
    errorbar((1:3)+dx,quart_pe_part(2,:),...
        quart_pe_part(2,:) - quart_pe_part(1,:),...
        quart_pe_part(3,:) - quart_pe_part(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3)-dx,quart_pe_pool(2,:),...
        quart_pe_pool(2,:) - quart_pe_pool(1,:),...
        quart_pe_pool(3,:) - quart_pe_pool(2,:),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    errorbar((1:3),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');

    plot([0,4],[pe0,pe0],'k:')

    title(['e_{PE} = ' num2str(dpe,'%0.1f') '\circ ; r_{PE} = ' num2str(cc.pe.r,'%0.2f')],...
      'FontSize',FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[27 52.9],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YMinorTick','on','FontSize',FontSize)

    subplot(122)
    errorbar((1:3)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3)-dx,quart_qe_pool(2,:),...
        quart_qe_pool(2,:) - quart_qe_pool(1,:),...
        quart_qe_pool(3,:) - quart_qe_pool(2,:),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    errorbar((1:3),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');
    plot([0,4],[qe0 qe0],'k:')

    title(['e_{QE} = ' num2str(dqe,'%0.1f') '% ; r_{QE} = ' num2str(cc.qe.r,'%0.2f')],...
      'FontSize',FontSize)
    ylabel('Quadrant Error (%)','FontSize',FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[0.1 49],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YAxisLocation','left',...
        'YMinorTick','on','FontSize',FontSize)

    l = legend('Part.','Pool','Actual');
    set(l,'Location','northwest','FontSize',FontSize-1)

    set(gcf,'PaperPosition',[1,1,10,3.5])
    
  end
end

%% ------ FIG 9 -----------------------------------------------------------
if flags.do_fig9
  
  % Model Settings
  latdivision = 0;            % lateral angle
  dlat = 10;

  % Experimental Settings
  Conditions = {'CL','N24','N9'};

  % Vocoder Settings 
  flow = 300;     % lowest corner frequency
  fhigh = 16000;  % highest corner frequency
  N = [inf,24,9];

  %% Computations
  s = data_baumgartner2014('pool');
  s = s(1); disp(['Listener: ' s.id]) % NH42: 8
  chance = [];
  for C = 1:length(Conditions)

    Cond = Conditions{C};

    %% Data

    % Experimental data
    data = data_goupell2010(Cond);
    for ll = 1:length(s)
      if sum(ismember({data.id},s(ll).id)) % if actual participant
        s(ll).mm1=data(ismember({data.id},s(ll).id)).mtx; 
        for ii = 1:length(latdivision)
          latresp = s(ll).mm1(:,7);
          idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
          mm2 = s(ll).mm1(idlat,:);
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

        if C==1
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


    %% Run Model

    for ll = 1:length(s)
      clear qe pe qe_t pe_t
      for ii = 1:length(latdivision)

        [p,rang] = baumgartner2014(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
              'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii});

        [ qe(ii),pe(ii) ] = pmv2ppp(p , s(ll).polang{ii} , rang);

        if sum(ismember({data.id},s(ll).id)) % if participant then actual targets
          [ qe_t(ii),pe_t(ii) ] = pmv2ppp( ...
              p , s(ll).polang{ii} , rang , s(ll).target{ii} );
        end

      end

      % Model results of pool
      s(ll).qe_pool(C,1) = mean(qe); 
      s(ll).pe_pool(C,1) = mean(pe);

      if sum(ismember({data.id},s(ll).id)) % if actual participant
        % Actual experimental results
        s(ll).qe_exp(C,1) = localizationerror(s(ll).mm1,'querrMiddlebrooks');
        s(ll).pe_exp(C,1) = localizationerror(s(ll).mm1,'rmsPmedianlocal');
        s(ll).Nt(C,1) = size(s(ll).mm1,1);
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

          if C==1; fp = figure; end
          subplot(1,3,C)
          plotbaumgartner2013(p,s(ll).polang{ii},rang,...
                s(ll).target{ii},s(ll).response{ii},...
                    'MarkerSize',MarkerSize,'cmax',0.05,'nocolorbar');
          title({['A: PE = ' num2str(s(ll).pe_exp(C,1),2) 'deg, QE = ' num2str(s(ll).qe_exp(C,1),2) '%'];...
            ['P: PE = ' num2str(s(ll).pe_part(C,1),2) 'deg, QE = ' num2str(s(ll).qe_part(C,1),2) '%']},'FontSize',FontSize-1)
          xlabel('Target Angle (deg)','FontSize',FontSize)
          ylabel('Response Angle (deg)','FontSize',FontSize)
          set(gca,'FontSize',FontSize-1)
          set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
          set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})

        end

      end

    end
  end
  
  varargout{1} = s;
  
end

%% ------ FIG 10 ----------------------------------------------------------
if flags.do_fig10
  
  fn = [mfilename('fullpath'),'_numchan.mat'];
  
  if amtredofile(fn,flags.redomode)
    
    autorefreshnotification(fn,flags)
    
    % Model Settings
    latdivision = 0; % lateral angle
    dlat = 10;

    % Experimental Settings
    Conditions = {'CL','N24','N18','N12','N9','N6','N3'};

    % Vocoder Settings 
    N = fliplr([3,6,9,12,18,24,30]);	% # of vocoder channels
    flow = 300;     % lowest corner frequency
    fhigh = 16000;  % highest corner frequency


    %% Computations
    s = data_baumgartner2014('pool');
    chance = [];
    for C = 1:length(Conditions)

      Cond = Conditions{C};

      %% Data

      % Experimental data
      data = data_goupell2010(Cond);
      for ll = 1:length(s)
        if sum(ismember({data.id},s(ll).id)) % if actual participant
          s(ll).mm1=data(ismember({data.id},s(ll).id)).mtx; 
          for ii = 1:length(latdivision)
            latresp = s(ll).mm1(:,7);
            idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
            mm2 = s(ll).mm1(idlat,:); 
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

          if C==1
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


      %% Run Model

      for ll = 1:length(s)
        clear qe pe qe_t pe_t
        for ii = 1:length(latdivision)

          [p,rang] = baumgartner2014(...
                s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
                'S',s(ll).S,'lat',latdivision(ii),...
                'polsamp',s(ll).polang{ii});

          [ qe(ii),pe(ii) ] = pmv2ppp(p , s(ll).polang{ii} , rang);

          if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
            [ qe_t(ii),pe_t(ii) ] = pmv2ppp( ...
                p , s(ll).polang{ii} , rang , s(ll).target{ii} );
          end

        end

        % Model results of pool
        s(ll).qe_pool(C,1) = mean(qe); 
        s(ll).pe_pool(C,1) = mean(pe);

        if sum(ismember({data.id},s(ll).id)) % if actual participant 
          % Actual experimental results
          s(ll).qe_exp(C,1) = localizationerror(s(ll).mm1,'querrMiddlebrooks');
          s(ll).pe_exp(C,1) = localizationerror(s(ll).mm1,'rmsPmedianlocal');
          s(ll).Nt(C,1) = size(s(ll).mm1,1);
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

    chance(:,8) = 240*rand(size(chance,1),1)-30;
    pe_chance = localizationerror(chance,'rmsPmedianlocal');
    qe_chance = localizationerror(chance,'querrMiddlebrooks');

    [r,p] =  corrcoef([s.qe_exp],[s.qe_part]);
    cc.qe.r = r(2);
    cc.qe.p = p(2);
    fprintf('QE: r = %0.2f, p = %0.3f\n',r(2),p(2));

    [r,p] =  corrcoef([s.pe_exp],[s.pe_part]);
    cc.pe.r = r(2);
    cc.pe.p = p(2);
    fprintf('PE: r = %0.2f, p = %0.3f\n',r(2),p(2));
    
    s = rmfield(s,{'spdtfs','spdtfs_c','Obj','mm1'});
    
    save(fn,'N', 'cc', 's', 'pe_chance', 'qe_chance',save_format);
  else
    load(fn);
  end
  varargout{1} = load(fn);
  
  %% Measures

  % Quartiles
  quart_pe_part = fliplr(quantile([s.pe_part]',[.25 .50 .75]));
  quart_qe_part = fliplr(quantile([s.qe_part]',[.25 .50 .75]));

  quart_pe_pool = fliplr(quantile([s.pe_pool]',[.25 .50 .75]));
  quart_qe_pool = fliplr(quantile([s.qe_pool]',[.25 .50 .75]));

  quart_pe_exp = fliplr(quantile([s.pe_exp]',[.25 .50 .75]));
  quart_qe_exp = fliplr(quantile([s.qe_exp]',[.25 .50 .75]));

  % RMS Differences
  % individual:
  Ntargets = [s.Nt]'; % # of targets
  relfreq = Ntargets/sum(Ntargets(:));
  sd_pe = ([s.pe_part]'-[s.pe_exp]').^2; % squared differences
  dpe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
  sd_qe = ([s.qe_part]'-[s.qe_exp]').^2;
  dqe = sqrt(relfreq(:)' * sd_qe(:));

  % Chance performance
  qe0 = qe_chance;
  pe0 = pe_chance;

  
  if flags.do_plot
    
    dx = 0.7;
    figure

    %% PE
    subplot(121)
    errorbar(fliplr(N)+dx,quart_pe_part(2,:),...
        quart_pe_part(2,:) - quart_pe_part(1,:),...
        quart_pe_part(3,:) - quart_pe_part(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar(fliplr(N)-dx,quart_pe_pool(2,:),...
        quart_pe_pool(2,:) - quart_pe_pool(1,:),...
        quart_pe_pool(3,:) - quart_pe_pool(2,:),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    errorbar(fliplr(N),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');
    plot([0,2*max(N)],[pe0,pe0],'k:')
    xlabel('Num. of Channels','FontSize',FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',FontSize)

    title(['e_{PE} = ' num2str(dpe,'%0.1f') '\circ ; r_{PE} = ' num2str(cc.pe.r,'%0.2f')],...
      'FontSize',FontSize)
    set(gca,'XLim',[1 32],'XTick',[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;9;12;18;24;'CL'},...
        'YLim',[29.1 54.5],...
        'YMinorTick','on','FontSize',FontSize)

    l = legend('Part.','Pool','Actual');
    set(l,'Location','southwest','FontSize',FontSize-2)

    %% QE
    subplot(122)
    plot([0,2*max(N)],[qe0,qe0],'k:')
    hold on
    errorbar(fliplr(N)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    errorbar(fliplr(N)-dx,quart_qe_pool(2,:),...
        quart_qe_pool(2,:) - quart_qe_pool(1,:),...
        quart_qe_pool(3,:) - quart_qe_pool(2,:),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    errorbar(fliplr(N),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');

    title(['e_{QE} = ' num2str(dqe,'%0.1f') '% ; r_{QE} = ' num2str(cc.qe.r,'%0.2f')],...
      'FontSize',FontSize)
    xlabel('Num. of Channels','FontSize',FontSize)
    ylabel('Quadrant Error (%)','FontSize',FontSize)
    set(gca,'XLim',[1 32],'XTick',[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;9;12;18;24;'CL'},...
        'YLim',[5.1 43],...
        'YMinorTick','on',...
        'YAxisLocation','left','FontSize',FontSize)

    set(gcf,'PaperPosition',[1,1,10,3.5])

  end
end

%% ------ FIG 11 ----------------------------------------------------------
if flags.do_fig11
  
  fn = [mfilename('fullpath'),'_nonindividual.mat'];
  
  if amtredofile(fn,flags.redomode)
    
    autorefreshnotification(fn,flags)
    
    % Settings
    latdivision = [-20,0,20];  % lateral center angles of SPs
    flow = 4e3;

    s = data_baumgartner2014('pool');
    ns = length(s);

    % DTFs of the SPs
    for ll = 1:ns
      for ii = 1:length(latdivision)

        s(ll).latang{ii} = latdivision(ii);
        s(ll).polangs{ii} = [];
        s(ll).spdtfs{ii} = [];
        [s(ll).spdtfs{ii},s(ll).polangs{ii}] = extractsp(...
            s(ll).latang{ii},s(ll).Obj);

      end
    end

    fprintf('\n Please wait a little! \n');
    qe = zeros(ns,ns,length(latdivision)); % init QEs
    pe = qe;           % init PEs
    pb = qe;           % init Polar Biases
    for ll = 1:ns    % listener
        for jj = 1:ns    % ears
            for ii = 1:length(latdivision) % SPs

              s(ll).p{jj,ii} = [];
              s(ll).respangs{ii} = [];
              [s(ll).p{jj,ii},s(ll).respangs{ii}] = baumgartner2014(...
                  s(jj).spdtfs{ii},s(ll).spdtfs{ii},s(ll).fs,...
                  'S',s(ll).S,'lat',s(ll).latang{ii},...
                  'polsamp',s(ll).polangs{ii},'flow',flow);

              [ qe(ll,jj,ii),pe(ll,jj,ii),pb(ll,jj,ii) ] = pmv2ppp( ...
                  s(ll).p{jj,ii} , s(jj).polangs{ii} , s(ll).respangs{ii});

            end
        end
        fprintf(' Subject %2u of %2u \n',ll,ns);
    end

    lat_weight = cos(pi*latdivision/180);     %lateral weight compensating compression of polar dimension
    lat_weight = lat_weight/sum(lat_weight);  % normalize
    lat_weight = repmat(reshape(lat_weight,[1,1,length(latdivision)]),[ns,ns,1]);
    qe_pool = sum(qe.*lat_weight,3);
    pe_pool = sum(pe.*lat_weight,3);
    pb_pool = sum(pb.*lat_weight,3);


    save(fn,'qe_pool', 'pe_pool','pb_pool',save_format);
  else
    load(fn);
  end
  varargout{1} = load(fn);
  
  data = data_middlebrooks1999;
  
  %% Model outcomes
  ns = size(pe_pool,1);
  own = eye(ns) == 1;
  other = not(own);
  pb_pool = abs(pb_pool);
  qe_own.quantiles = quantile(qe_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pe_own.quantiles = quantile(pe_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pb_own.quantiles = quantile(pb_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
  qe_own.mean = mean(qe_pool(own));
  pe_own.mean = mean(pe_pool(own));
  pb_own.mean = mean(pb_pool(own));
  
  qe_other.quantiles = quantile(qe_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pe_other.quantiles = quantile(pe_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pb_other.quantiles = quantile(pb_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
  qe_other.mean = mean(qe_pool(other));
  pe_other.mean = mean(pe_pool(other));
  pb_other.mean = mean(pb_pool(other));
  
  if flags.do_plot
    dx = -0.2;
    Marker = 'ks';
    data.Marker = 'ko';
    MFC = 'k'; % Marker Face Color
    data.MFC = 'w';
    
    fig = figure;
    subplot(131)
    middlebroxplot(1-dx,qe_own.quantiles,MarkerSize)
    plot(1-dx,qe_own.mean,Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',MFC)
    middlebroxplot(1+dx,data.qe_own.quantiles,MarkerSize)
    plot(1+dx,data.qe_own.mean,data.Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',data.MFC)
    middlebroxplot(2-dx,qe_other.quantiles,MarkerSize)
    plot(2-dx,qe_other.mean,Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',MFC)
    middlebroxplot(2+dx,data.qe_other.quantiles,MarkerSize)
    plot(2+dx,data.qe_other.mean,data.Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',data.MFC)
    ylabel('Quadrant Errors (%)','FontSize',FontSize)
    set(gca,'YLim',[-2 43],'XLim',[0.5 2.5],...
      'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',FontSize)

    subplot(132)
    plot(1-dx,pe_own.mean,Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',MFC)
    hold on
    plot(1+dx,data.pe_own.mean,data.Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',data.MFC)

    middlebroxplot(1-dx,pe_own.quantiles,MarkerSize)
    plot(1-dx,pe_own.mean,Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',MFC)
    middlebroxplot(1+dx,data.pe_own.quantiles,MarkerSize)
    plot(1+dx,data.pe_own.mean,data.Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',data.MFC)
    middlebroxplot(2-dx,pe_other.quantiles,MarkerSize)
    plot(2-dx,pe_other.mean,Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',MFC)
    middlebroxplot(2+dx,data.pe_other.quantiles,MarkerSize)
    plot(2+dx,data.pe_other.mean,data.Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',data.MFC)
    ylabel('Local Polar RMS Error (deg)','FontSize',FontSize)
    set(gca,'YLim',[-2 62],'XLim',[0.5 2.5],...
      'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',FontSize)

    subplot(133)
    middlebroxplot(1-dx,pb_own.quantiles,MarkerSize)
    plot(1-dx,pb_own.mean,Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',MFC)
    middlebroxplot(1+dx,data.pb_own.quantiles,MarkerSize)
    plot(1+dx,data.pb_own.mean,data.Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',data.MFC)
    middlebroxplot(2-dx,pb_other.quantiles,MarkerSize)
    plot(2-dx,pb_other.mean,Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',MFC)
    middlebroxplot(2+dx,data.pb_other.quantiles,MarkerSize)
    plot(2+dx,data.pb_other.mean,data.Marker,'MarkerSize',MarkerSize,'MarkerFaceColor',data.MFC)

    ylabel('Magnitude of Elevation Bias (deg)','FontSize',FontSize)
    set(gca,'YLim',[-2 55],'XLim',[0.5 2.5],...
      'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',FontSize)
  end
end

%% ------ FIG 12 ----------------------------------------------------------
if flags.do_fig12
  
  fn = [mfilename('fullpath'),'_ripples.mat'];
  
  if amtredofile(fn,flags.redomode)
    
    autorefreshnotification(fn,flags)
    
    do_exp1 = true;
    do_exp2 = true;
    plotpmv = false;

    density = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8]; % ripples/oct
    depth =   10:10:40;        % ripple depth (peak-to-trough) in dB

    %% Stimulus: 
    % 250-ms bursts, 20-ms raised-cosine fade in/out, flat from 0.6-16kHz

    fs = 48e3;    % sampling rate
    flow = 1e3;   % lower corner frequency of ripple modification in Hz
    fhigh = 16e3; % upper corner frequency of ripple modification in Hz
    Nf = 2^10;    % # Frequency bins

    f = 0:fs/2/Nf:fs/2;	% frequency bins
    id600 = find(f<=600,1,'last'); % index of 600 Hz (lower corner frequency of stimulus energy)
    idlow = find(f<=flow,1,'last'); % index of flow (ripples)
    idhigh = find(f>=fhigh,1,'first');  % index of fhigh (ripples)
    N600low = idlow - id600 +1;   % # bins without ripple modification
    Nlowhigh = idhigh - idlow +1; % # bins with ripple modification     % 
    O = log2(f(idlow:idhigh)/1e3);   % freq. trafo. to achieve equal ripple density in log. freq. scale

    % Raised-cosine "(i.e., cos^2)" ramp 1/8 octave wide
    fup = f(idlow)*2^(1/8);       % upper corner frequency of ramp upwards 
    idup = find(f<=fup,1,'last');
    Nup = idup-idlow+1;
    rampup = cos(-pi/2:pi/2/(Nup-1):0).^2;
    fdown = f(idhigh)*2^(-1/8);  % lower corner frequency of ramp downwards
    iddown = find(f>=fdown,1,'first');
    Ndown = idhigh-iddown+1;
    rampdown = cos(0:pi/2/(Ndown-1):pi/2).^2;
    ramp = [rampup ones(1,Nlowhigh-Nup-Ndown) rampdown];
    ramp = [-inf*ones(1,id600-1) zeros(1,N600low) ramp -inf*ones(1,Nf - idhigh)];

    % Ripples of Experiment I
    Sexp1 = zeros(Nf+1,length(density),2);  % 3rd dim: 1:0-phase 2:pi-phase
    Sexp1(idlow:idhigh,:,1) = (40/2* sin(2*pi*density'*O+ 0))';  % depth: 40dB, 0-phase
    Sexp1(idlow:idhigh,:,2) = (40/2* sin(2*pi*density'*O+pi))';  % depth: 40dB, pi-phase
    Sexp1 = repmat(ramp',[1,length(density),2]) .* Sexp1;
    Sexp1 = [Sexp1;Sexp1(Nf-1:-1:2,:,:)];
    Sexp1(isnan(Sexp1)) = -100;
    sexp1 = ifftreal(10.^(Sexp1/20),2*Nf);
    sexp1 = circshift(sexp1,Nf);  % IR corresponding to ripple modification
    % Ripples of Experiment II
    Sexp2 = zeros(Nf+1,length(depth),2);  % 3rd dim: 1:0-phase 2:pi-phase
    Sexp2(idlow:idhigh,:,1) = (depth(:)/2*sin(2*pi*1*O+ 0))';  % density: 1 ripple/oct, 0-phase
    Sexp2(idlow:idhigh,:,2) = (depth(:)/2*sin(2*pi*1*O+pi))';  % density: 1 ripple/oct, pi-phase
    Sexp2 = repmat(ramp',[1,length(depth),2]) .* Sexp2;
    Sexp2 = [Sexp2;Sexp2(Nf-1:-1:2,:,:)];
    Sexp2(isnan(Sexp2)) = -100;
    sexp2 = ifftreal(10.^(Sexp2/20),2*Nf);
    sexp2 = circshift(sexp2,Nf);  % IR corresponding to ripple modification


    %% Modeling
    for do = 0:1

      if do == 1
        s = data_baumgartner2014('pool');
      else % recalib
        s = data_baumgartner2014('pool','recalib','do',0);
      end

    latseg = 0;   % centers of lateral segments
    runs = 5;     % # runs of virtual experiments

    pe_exp1 = zeros(length(latseg),length(s),length(density),2);
    pe_exp2 = zeros(length(latseg),length(s),length(depth),2);
    pe_flat = zeros(length(latseg),length(s));
    for ss = 1:length(s)
      for ll = 1:length(latseg)

        [spdtfs,polang] = extractsp(latseg(ll),s(ss).Obj);

        % target elevation range of +-60°
        idt = find( polang<=60 | polang>=120 );
        targets = spdtfs(:,idt,:);
        tang = polang(idt);

        [pflat,rang] = baumgartner2014(targets,spdtfs,...
            'S',s(ss).S,'polsamp',polang,...
            'lat',latseg(ll),'stim',[1;0],'do',do); % Impulse
        mflat = virtualexp(pflat,tang,rang,'runs',runs);
        pe_flat(ll,ss) = pemacpherson2003(mflat,mflat);

        if plotpmv; plotbaumgartner2013(pflat,tang,rang,mflat(:,6),mflat(:,8));title(num2str(pe_flat(ll,ss),2));pause(0.5); end 

        if do_exp1  % Exp. I
        for ii = 1:2*length(density)

          [p,rang] = baumgartner2014(targets,spdtfs,...
            'S',s(ss).S,'polsamp',polang,...
            'lat',latseg(ll),'stim',sexp1(:,ii),'do',do);
          m = virtualexp(p,tang,rang,'runs',runs);
          pe_exp1(ll,ss,ii) = pemacpherson2003(m,mflat);% - pe_flat(ll,ss);

          if plotpmv; plotbaumgartner2013(p,tang,rang,m(:,6),m(:,8));title([num2str(density(mod(ii-1,10)+1)) 'ripples/oct; PE:' num2str(pe_exp1(ll,ss,ii),2) '%']);pause(0.5); end

        end
        end

        if do_exp2 % Exp. II
        for ii = 1:2*length(depth)

          [p,rang] = baumgartner2014(targets,spdtfs,...
            'S',s(ss).S,'polsamp',polang,...
            'lat',latseg(ll),'stim',sexp2(:,ii),'do',do);
          m = virtualexp(p,tang,rang,'runs',runs);
          pe_exp2(ll,ss,ii) = pemacpherson2003(m,mflat);% - pe_flat(ll,ss);

          if plotpmv; plotbaumgartner2013(p,tang,rang,m(:,6),m(:,8));title([num2str(depth(mod(ii-1,4)+1)) 'dB; PE:' num2str(pe_exp2(ll,ss,ii),2) '%']);pause(0.5); end

        end
        end

      end

      fprintf('\n %2.0f of %2.0f subjects completed',ss,length(s))

    end
    fprintf('\n')

    if length(latseg) > 1
      pe_exp1 = squeeze(mean(pe_exp1));
      pe_exp2 = squeeze(mean(pe_exp2));
      pe_flat = squeeze(mean(pe_flat));
    else 
      pe_exp1 = squeeze(pe_exp1);
      pe_exp2 = squeeze(pe_exp2);
      pe_flat = squeeze(pe_flat);
    end

    %% Save
      if do==0
        noDCN.pe_exp1 = pe_exp1;
        noDCN.pe_exp2 = pe_exp2;
        noDCN.pe_flat = pe_flat;
        delete(which('baumgartner2014calibration.mat'))
      end
    end

    save(fn,'pe_exp1','pe_exp2','pe_flat','noDCN',save_format);
  else
    load(fn);
  end
  varargout{1} = load(fn);
  
  dcn_flag = true;
  
  % Original data:
  data = data_macpherson2003;


  %% Phase condition handling
  pe_exp1 = mean(pe_exp1,3);
  data.pe_exp1 = mean(data.pe_exp1,3);
  pe_exp2 = mean(pe_exp2,3);
  data.pe_exp2 = mean(data.pe_exp2,3);
  if dcn_flag
      noDCN.pe_exp1 = mean(noDCN.pe_exp1,3);
      noDCN.pe_exp2 = mean(noDCN.pe_exp2,3);
  end
  idphase = 1;


  %% Statistics
  quart_pe_flat = quantile(pe_flat,[.25 .50 .75]);
  quart_pe_data_flat = quantile(data.pe_flat,[.25 .50 .75]);

  quart_pe_exp1 = quantile(pe_exp1,[.25 .50 .75]);
  quart_pe_data_exp1 = quantile(data.pe_exp1,[.25 .50 .75]);

  quart_pe_exp2 = quantile(pe_exp2,[.25 .50 .75]);
  quart_pe_data_exp2 = quantile(data.pe_exp2,[.25 .50 .75]);

  if dcn_flag
      noDCN.quart_pe_flat = quantile(noDCN.pe_flat,[.25 .50 .75]);
      noDCN.quart_pe_exp1 = quantile(noDCN.pe_exp1,[.25 .50 .75]);
      noDCN.quart_pe_exp2 = quantile(noDCN.pe_exp2,[.25 .50 .75]);
  end

  
  if flags.do_plot
    
    dx = 1.05;
    % Exp1
    fig = figure;
    subplot(5,1,1:3)
    errorbar(data.density*dx,quart_pe_exp1(2,:,idphase),...
        quart_pe_exp1(2,:,idphase) - quart_pe_exp1(1,:,idphase),...
        quart_pe_exp1(3,:,idphase) - quart_pe_exp1(2,:,idphase),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    if dcn_flag
        errorbar(data.density,noDCN.quart_pe_exp1(2,:,idphase),...
        noDCN.quart_pe_exp1(2,:,idphase) - noDCN.quart_pe_exp1(1,:,idphase),...
        noDCN.quart_pe_exp1(3,:,idphase) - noDCN.quart_pe_exp1(2,:,idphase),...
        'kd--','MarkerSize',MarkerSize-1,...
        'MarkerFaceColor','k');
    end
    errorbar(data.density/dx,quart_pe_data_exp1(2,:,idphase),...
        quart_pe_data_exp1(2,:,idphase) - quart_pe_data_exp1(1,:,idphase),...
        quart_pe_data_exp1(3,:,idphase) - quart_pe_data_exp1(2,:,idphase),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');
    set(gca,'XScale','log','YMinorTick','on')
    set(gca,'XLim',[0.25/1.2 8*1.2],'XTick',data.density,'YLim',[-1 75],'FontSize',FontSize)
    xlabel('Ripple Density (ripples/octave)','FontSize',FontSize)
    ylabel({'Polar Error Rate (%)'},'FontSize',FontSize)

    if dcn_flag
        leg = legend('P (with DCN)','P (w/o DCN)','A');
    else
        leg = legend('Predicted','Actual');
    end
    set(leg,'FontSize',FontSize-2,'Location','northeast')

    %% Exp2
    data.depth = [0 data.depth];
    quart_pe_exp2 = [quart_pe_flat(:) , quart_pe_exp2];
    if dcn_flag
      noDCN.quart_pe_exp2 = [noDCN.quart_pe_flat(:) , noDCN.quart_pe_exp2];
    end
    quart_pe_data_exp2 = [quart_pe_data_flat(:) , quart_pe_data_exp2];
    
    subplot(5,1,4:5)
    errorbar(data.depth-1,quart_pe_exp2(2,:,idphase),...
        quart_pe_exp2(2,:,idphase) - quart_pe_exp2(1,:,idphase),...
        quart_pe_exp2(3,:,idphase) - quart_pe_exp2(2,:,idphase),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    if dcn_flag
        errorbar(data.depth,noDCN.quart_pe_exp2(2,:,idphase),...
        noDCN.quart_pe_exp2(2,:,idphase) - noDCN.quart_pe_exp2(1,:,idphase),...
        noDCN.quart_pe_exp2(3,:,idphase) - noDCN.quart_pe_exp2(2,:,idphase),...
        'kd--','MarkerSize',MarkerSize-1,...
        'MarkerFaceColor','k');
    end
    errorbar(data.depth+1,quart_pe_data_exp2(2,:,idphase),...
        quart_pe_data_exp2(2,:,idphase) - quart_pe_data_exp2(1,:,idphase),...
        quart_pe_data_exp2(3,:,idphase) - quart_pe_data_exp2(2,:,idphase),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');
    set(gca,'XLim',[data.depth(1)-5 data.depth(end)+5],'XTick',data.depth,'YLim',[-4 69],'YMinorTick','on','FontSize',FontSize)
    xlabel('Ripple Depth (dB)','FontSize',FontSize)
    ylabel({'Polar Error Rate (%)'},'FontSize',FontSize)
    
  end
end

%% ------ FIG 13 ----------------------------------------------------------
if flags.do_fig13
  
  fn = [mfilename('fullpath'),'_highfreqatten.mat'];
  
  if amtredofile(fn,flags.redomode)
    
    if not(exist('HarvardWords','dir'))
      disp('The Harvard word list is missing.') 
      disp('Please, contact Virginia Best (ginbest@bu.edu) or Craig Jin (craig.jin@sydney.edu.au) for providing their speech recordings.')
      disp(['Then, move the folder labeled HarvardWords to: ' fullfile(amtbasepath,'signals') '.'])
      return
    end
    
    autorefreshnotification(fn,flags)
    
    tempfn = fullfile(amtbasepath,'experiments','exp_baumgartner2014_highfreqatten'); % temporary folder
    mkdir(tempfn)
    
    %% Settings
    latseg = 0;%[-20,0,20];   % centers of lateral segments
    NsampModel = 260; % # of modeled speech samples (takes 30min/sample)
    startSamp = 1;

    saveflag = true;
    plotpmv = false;
    plotspec = false;


    %% Loda Data

    % Listeners
    s = data_baumgartner2014('pool');

    % Speech Samples from Harvard Word list
    fs_orig = 80e3; % Hz
    fs = s(1).fs;   % Hz
    p_resamp = fs/fs_orig;
    kk = 1;
    if NsampModel <= 51
      Nsamp = NsampModel;
      Nlists = 1;
    else
      Nsamp = 260;
      Nlists = 5;
    end
    lsamp = 120000*p_resamp;
    speechsample = cell(Nsamp,1);
    for ii = 1:Nlists
      tmp.path = fullfile(mfilename,'HarvardWords',['list' num2str(ii,'%1.0u')]);
      tmp.dir = dir(fullfile(tmp.path,'*.mat'));
      for jj = 1:length(tmp.dir)
        if jj > Nsamp; break; end
        load(fullfile(tmp.path,tmp.dir(jj).name))
        signal = resample(word,p_resamp*10,10);
        env = filter(normpdf(0:0.001:10,0,1),1,signal.^2);
        idon = max(find(env > 5e7,1,'first')-1e3,1);
        idoff = min(find(env > 5e7,1,'last')+1e3,lsamp);
        lwin = idoff-idon+1;
        speechsample{kk} = signal(idon:idoff) .* tukeywin(lwin,0.01)';
        kk = kk + 1;
      end
    end


    % FIR Low-pass filters at 8kHz
    % Brick-wall (aka sinc-filter): fir1(200,1/3) -> -60 dB
    fn_filters = 'exp_baumgartner2014_highfreqatten_filters.mat';
    if not(exist(fn_filters,'file'))
      disp(['Downloading ' fn_filters ' from http://www.kfs.oeaw.ac.at/']);
      targetfn = fullfile(amtbasepath,'humandata',fn_filters);
      sourcefn = ['http://www.kfs.oeaw.ac.at/research/experimental_audiology/projects/amt/' fn_filters];
      urlwrite(sourcefn,targetfn);
    end
    load(fn_filters)
    lp{1} = [1 zeros(1,100)];
    lp{2} = fir20db;
    lp{3} = fir40db;
    lp{4} = fir60db;

    %% Model Data

    ape_BBnoise = zeros(1,length(s),length(latseg));
    qe_BBnoise = ape_BBnoise;
    for ss = 1:length(s)
      for ll = 1:length(latseg)
        [spdtfs,polang] = extractsp(latseg(ll),s(ss).Obj);
        [p,rang] = baumgartner2014(spdtfs,spdtfs,...
              'S',s(ss).S,'polsamp',polang,'lat',latseg(ll),'notprint');
        ape_BBnoise(1,ss,ll) = apebest2005(p,polang,rang);
        qe_BBnoise(1,ss,ll) = pmv2ppp(p,polang,rang);

        if plotpmv; figure; plotbaumgartner2013(p,polang,rang); title(num2str(ape_BBnoise(1,ss,ll),2)); end

      end
    end

    %%
    ape_all = zeros(length(lp),length(s),length(latseg));
    qe_all = ape_all;
    for kk = startSamp:NsampModel % start with 104
      for ss = 1:length(s)
        for ll = 1:length(latseg)
          for ii = 1:length(lp)

            stim = filter(lp{ii},1,speechsample{kk});

            if plotspec; figure; audspecgram(stim(:),fs,'dynrange',150); end

            [spdtfs,polang] = extractsp(latseg(ll),s(ss).Obj);
            [p,rang] = baumgartner2014(spdtfs,spdtfs,...
              'S',s(ss).S,'polsamp',polang,...
              'lat',latseg(ll),'stim',stim,'notprint');
            ape_all(ii,ss,ll) = apebest2005(p,polang,rang);
            qe_all(ii,ss,ll) = pmv2ppp(p,polang,rang);

            if plotpmv; figure; plotbaumgartner2013(p,polang,rang); title(num2str(ape_all(ii,ss,ll),2)); end

          end
        end
      end

      % % Pool Lateral Segments
      if length(latseg) > 1
        ape_BBnoise = mean(ape_BBnoise,3);
        ape_all = mean(ape_all,3);
        qe_BBnoise = mean(qe_BBnoise,3);
        qe_all = mean(qe_all,3);
      end

      if saveflag
        savename = fullfile(tempfn,['result_best2005speech_samp' num2str(kk)]);
        save(savename,'ape_all','qe_all','ape_BBnoise','qe_BBnoise')
      end
      fprintf('%1.0u of %2.0u samples completed\n',kk,NsampModel)
    end
    
    results = dir(fullfile(tempfn,'result_best2005speech_samp*.mat'));
    tmp = load(fullfile(tempfn,results(1).name));
    ape_BBnoise = tmp.ape_BBnoise;
    qe_BBnoise = tmp.qe_BBnoise;
    ape_all = zeros([size(tmp.ape_all) length(results)]);
    qe_all = ape_all;
    for ii = 1:length(results)
      tmp = load(fullfile(tempfn,results(ii).name));
      ape_all(:,:,ii) = tmp.ape_all;
      qe_all(:,:,ii) = tmp.qe_all;
    end
    
    save(fn,'ape_all','qe_all','ape_BBnoise','qe_BBnoise',save_format);
%     rmdir(tempfn,'s');
    
  else
    load(fn);
  end
  varargout{1} = load(fn);
  
  data = data_best2005;

  % Pool Samples
  ape_pooled = mean(ape_all,3);
  qe_pooled = mean(qe_all,3);

  % Confidence Intervals or standard errors
  df_speech = size(ape_all,2)-1;%*size(ape_all,3)-1;
  tquant_speech = 1;%icdf('t',.975,df_speech);
  seape_speech = std(ape_pooled,0,2)*tquant_speech/(df_speech+1);
  df_noise = size(ape_BBnoise,2)-1;
  tquant_noise = 1;%icdf('t',.975,df_noise);
  seape_noise = std(ape_BBnoise,0,2)*tquant_noise/(df_noise+1);
  seape = [seape_noise;seape_speech];

  % Means
  ape = mean([ape_BBnoise ; ape_pooled],2);
  qe = mean([qe_BBnoise ; qe_pooled],2);


  if flags.do_plot
    
    dx = 0;
    xticks = 0:size(ape_all,1);
    fig = figure;
    subplot(211)
    h(1) = errorbar(xticks-dx,ape,seape,'ko');
    set(h(1),'MarkerFaceColor','k','MarkerSize',MarkerSize,'LineStyle','-')
    hold on
    h(2) = errorbar(xticks+dx,data.ape,data.seape,'ko');
    set(h(2),'MarkerFaceColor','w','MarkerSize',MarkerSize,'LineStyle','-')
    ylabel('| \theta - \vartheta | (deg)','FontSize',FontSize)
    set(gca,'XTick',xticks,'XTickLabel',[],'FontSize',FontSize)
    set(gca,'XLim',[-0.5 4.5],'YLim',[14 79],'YMinorTick','on')

    pos = get(gca,'Position');
    pos(2) = pos(2)-0.1;
    set(gca,'Position',pos)

    subplot(212)
    h(1) = plot(xticks-dx,qe,'ko');
    set(h(1),'MarkerFaceColor','k','MarkerSize',MarkerSize,'LineStyle','-')
    hold on
    h(2) = plot(xticks([1 2 5])+dx,data.qe([1 2 5]),'ko');
    set(h(2),'MarkerFaceColor','w','MarkerSize',MarkerSize,'LineStyle','-')
    ylabel('QE (%)','FontSize',FontSize)
    set(gca,'XTick',xticks,'XTickLabel',data.meta,'FontSize',FontSize,...
      'XLim',[-0.5 4.5],'YLim',[-5 45],'YMinorTick','on')

  end
end
  
end



%% ------------------------------------------------------------------------
%  ---- INTERNAL FUNCTIONS ------------------------------------------------
%  ------------------------------------------------------------------------
function autorefreshnotification(fn,flags)
if flags.do_autorefresh
  disp(['Calculation of ' fn ' started.'])
  disp('Results can be also downloaded here:') 
  disp(' http://www.kfs.oeaw.ac.at/research/experimental_audiology/projects/amt/exp_baumgartner2014.zip')
  disp('Unzip the folder and move the single files into:')
  disp([' ' fullfile(amtbasepath,'experiments')])
end
end

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
        hwinfade = fade(hwin,512,24,96,192);
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
N = duration*GaussRate; % number of pulses
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

function middlebroxplot(x,quantiles,MarkerSize)

lilen = 0.14; % length of horizontal lines

% Symbols
plot(x,quantiles(1),'kx','MarkerSize',MarkerSize) % min
hold on
plot(x,quantiles(7),'kx','MarkerSize',MarkerSize) % max

% Horizontal lines
line(x+0.5*[-lilen,lilen],repmat(quantiles(2),2),'Color','k') % lower whisker
line(x+[-lilen,lilen],repmat(quantiles(3),2),'Color','k') % 25% Quartile
line(x+[-lilen,lilen],repmat(quantiles(4),2),'Color','k') % Median
line(x+[-lilen,lilen],repmat(quantiles(5),2),'Color','k') % 75% Quartile
line(x+0.5*[-lilen,lilen],repmat(quantiles(6),2),'Color','k') % upper whisker

% Vertical lines
line([x,x],quantiles(2:3),'Color','k') % connector lower whisker
line([x,x],quantiles(5:6),'Color','k') % connector upper whisker
line([x,x]-lilen,quantiles([3,5]),'Color','k') % left box edge
line([x,x]+lilen,quantiles([3,5]),'Color','k') % left box edge

end

function per = pemacpherson2003(mmod,mflat)
%PEMACPHERSON2003 polar error rate used in Macpherson & Middlebrooks (2000)
%   Usage: per = pemacpherson2003(mmod,mflat)
%
%   Input parameters:
%     mmod  : localization responses to modified spectra
%     mflat : localization responses to flat spectra (baseline)
%
%   Output parameters:
%     per   : polar error rate in %
%
%   `pemacpherson2003()` assesses localization accuracy for flat- and 
%   shaped-spectrum targets by measuring the deviation of responses from 
%   the linear predictors obtained by and ad hoc selective, iterative 
%   regression procedure. Polar errors are defined by showing a deviation 
%   of >45° with respect to the linear flat stimulus prediction.
%
%   See also: sirproceduremacpherson2000

% AUTHOR: Robert Baumgartner

plotflag = false;
tol = 45; % tolerance of polar angle to be not counted as polar error

mflatf = mflat( abs(mflat(:,7))<=30 & mflat(:,6)< 90 ,:); % frontal central data only
mflatr = mflat( abs(mflat(:,7))<=30 & mflat(:,6)>=90 ,:); % rear central data only

mmodf = mmod( abs(mmod(:,7))<=30 & mmod(:,6)< 90 ,:); % frontal central data only
mmodr = mmod( abs(mmod(:,7))<=30 & mmod(:,6)>=90 ,:); % rear central data only

[f,r] = sirproceduremacpherson2000(mflat);
yf=f.b(2)*mflatf(:,6)+f.b(1); % linear prediction for flat, frontal targets
yr=r.b(2)*mflatr(:,6)+r.b(1); % linear prediction for flat, rear targets

f.pe = sum(abs(wrapTo180(mmodf(:,8)-yf)) > tol);
r.pe = sum(abs(wrapTo180(mmodr(:,8)-yr)) > tol);

per = (f.pe + r.pe) / size(mmod,1) * 100; % polar error rate

if plotflag 
  
  Loca_PlotResponse(mmod(:,6),mmod(:,8),'polar');
  hold on
  plot(mflatf(:,6),yf,'LineWidth',2);
  plot(mflatr(:,6),yr,'LineWidth',2);
  plot(mflatf(:,6),yf+45,'r','LineWidth',2);
  plot(mflatf(:,6),yf-45,'r','LineWidth',2);
  plot(mflatr(:,6),yr+45,'r','LineWidth',2);
  plot(mflatr(:,6),yr-45,'r','LineWidth',2);
  title(['PER: ' num2str(per,2) '%'])
  
end

end

function [f,r] = sirproceduremacpherson2000(m)
%SIRPROCEDUREMACPHERSON2003 ad hoc selective, iterative regression
%procedure proposed in Macpherson & Middlebrooks (2000)
%   Usage: [f,r] = sirproceduremacpherson2000(m)
%
%   Input parameters:
%     m: data matrix
%
%   Output parameters:
%     f : structure containing regression results the procedure converged 
%         to for the frontal hemisphere. For detailed description of the
%         fields see help of `regress()`
%     r : structure for rear hemisphere
%
%   `sirproceduremacpherson2000()` performs an ad hoc selective, iterative
%   regression procedure in order to exclude outliers and reversals and 
%   isolate the main concentration of responses in the computation of the 
%   linear fits. Outlier distance criterion: 40°
%
%   See also: regress

% AUTHOR: Robert Baumgartner

delta = 40; % outlier tolerance in degrees
plotflag = false;

%% Front
mf = m( abs(m(:,7))<=30 & m(:,6)<90 ,:); % frontal central data only
idf = find(mf(:,8)<90);   % indices correct forntal responses (init)
if plotflag; Loca_PlotResponse(mf(:,6),mf(:,8),'polar'); end

if length(idf)<2
  f.b = zeros(1,2);
  f.stats = zeros(1,4);
else
  old=[];
  while not(isequal(idf,old))
    [f.b,f.bint,f.r,f.rint,f.stats]=regress(mf(idf,8),[ones(length(idf),1) mf(idf,6)]);
    y=f.b(2)*mf(:,6)+f.b(1);
    if plotflag;  hold on; plot(mf(:,6),y); end
    old = idf;
    idf=find( abs(wrapTo180(mf(:,8)-y)) < delta); 
  end	
end

%% Rear
mr = m( abs(m(:,7))<=30 & m(:,6)>=90 ,:); % rear central data only
idr = find(mr(:,8)>=90);  % indices correct rear responses (init)
if plotflag; Loca_PlotResponse(mr(:,6),mr(:,8),'polar'); end

if length(idr)<2
  r.b = zeros(1,2);
  r.stats = zeros(1,4);
else
  old=[];
  while not(isequal(idr,old))	
    [r.b,r.bint,r.r,r.rint,r.stats]=regress(mr(idr,8),[ones(length(idr),1) mr(idr,6)]);
    y=r.b(2)*mr(:,6)+r.b(1);
    if plotflag;  hold on; plot(mr(:,6),y); end
    old = idr;
    idr=find( abs(wrapTo180(mr(:,8)-y)) < delta);
  end	
end

end

function m = virtualexp(p,tang,rang,varargin)
%VIRTUALEXP Response patterns of virtual localization experiments
%   Usage:    m = virtualexp(p,tang,rang)
%
%   Input parameters:
%     p       : prediction matrix containing probability mass vectors (PMVs) 
%               for the polar response angle as a function of the polar  
%               target angle (1st dim: response angle, 2nd dim: target
%               angle)
%     rang    : polar response angles
%     tang    : polar target angles
%
%   Output parameter:
%     m       : item list of virtual experiment
%               Columns: 
%                1:4 ...   azi_target,ele_target,azi_response,ele_response
%                5:8 ...   lat_target,pol_target,lat_response,pol_response
%                9   ...   F/B-C resolved pol_response
%
%   `virtualexp(...)` runs virtual localization experiments where the
%   response behavior is based on (predicted) polar response PMVs.
%
%   `virtualexp` accepts the following optional parameters:
%
%     'runs',runs    	Define the number of runs. 
%                    	Default value is 1.
%
%     'targetset',ts  Define the set of polar target angles.
%                    	As default 'tang' is used.
%
%     'lat',lat     	Define the lateral target angles. 
%                    	Default value is 0°.
%
%   See also: baumgartner2013, plotbaumgartner2013
%
%   References: baumgartner2013assessment baumgartner2012modelling langendijk2002contribution patterson1988efficient dau1996qmeI

    
% AUTHOR: Robert Baumgartner

definput.keyvals.runs = 1;
definput.keyvals.targetset = [];
definput.keyvals.lat = 0;
% definput.flags.colorbar = {'colorbar','nocolorbar'};
[flags,kv]=ltfatarghelper({'runs','targetset'},definput,varargin);

if isempty(kv.targetset)
  kv.targetset = tang;
end


%% Run experiments
nt=length(kv.targetset);
m = nan(nt*kv.runs,9);
m(:,5) = kv.lat;
m(:,6) = repmat(kv.targetset(:),kv.runs,1);
m(:,7) = kv.lat;
% kv.targetset = round(kv.targetset);
if length(tang) > 1
  tangbound = tang(:)+0.5*diff([tang(1)-diff(tang(1:2));tang(:)]);
else
  tangbound = tang;
end
post=zeros(nt,1); % indices of target positions
for ii = 1:nt
    post(ii) = find(tangbound>=kv.targetset(ii),1);
end

posr=zeros(nt,1);
for rr=1:kv.runs
  for jj = 1:nt 
    posr(jj) = discreteinvrnd(p(:,post(jj)),1);
    m(jj+(rr-1)*nt,8) = rang(posr(jj));
  end
end



end

function [ X ] = discreteinvrnd(p,n,m)
% DISCRETEINVRND implements an inversion method for a discrete distribution
% with probability mass vector p for n trials
% Usage:    [ X ] = discreteinvrnd(p,n)
%
% AUTHOR : Robert Baumgartner

if ~exist('m','var')
    m=1;
end

p = p/sum(p);   % ensure probability mass vector
c = cumsum(p);
t = max(c)*rand(n,m); % rand returns ]0,1]
X = zeros(n,m);
for jj = 1:m
    for ii = 1:n
        X(ii,jj) = find(c >= t(ii,jj) ,1);
    end
end

end

function ape = apebest2005(p,tang,rang)
%APEBEST2005 PMV to PPP conversion
%   Usage:  [ ape ] = apebest2005( p,tang,rang );
%
%   Input parameters:
%     p          : prediction matrix (response PMVs)
%     tang       : possible polar target angles. As default, ARI's MSP 
%                  polar angles in the median SP is used.
%     rang       : polar angles of possible response angles.
%                  As default regular 5°-sampling is used (-30:5:210).    
%
%   Output parameter:
%     ape        : absolute polar angle error in degrees
%

% AUTHOR : Robert Baumgartner

nt = length(tang);

apet = zeros(nt,1);
for ii = 1:nt % for all target positions
  
    d = tang(ii)-rang;                 % wraped angular distance between tang & rang
    iduw = (d < -180) | (180 < d);     % 180°-unwrap indices
    d(iduw) = mod(d(iduw) + 180,360) - 180; % 180° unwrap
    d = abs(d);                        % absolut distance
    apet(ii) = sum( p(:,ii) .* d');     % absolut polar angle error for target ii
    
end

ape = mean(apet);

end
