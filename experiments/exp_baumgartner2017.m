function data = exp_baumgartner2017(varargin)
%EXP_BAUMGARTNER2017 - Experiments of Baumgartner et al. (2017).
%   Usage: data = exp_baumgartner2017(flag) 
%
%   `exp_baumgartner2017(flag)` reproduces figures of the study from 
%   Baumgartner et al. (2016).
%
% %   Optional fields of output *data* structure:
% %
% %   `data.contralateralGain`
% %      contralateral gain of binaural weighting function
%
%
%   The following flags can be specified
%
%     `boyd2012`    
%         Model experiments from Boyd et al. (2012; Fig.1, top):
%         Average externalization ratings of 1 talker for NH participants 
%         against mix point as a function of microphone position (ITE/BTE) 
%         and frequency response (BB/LP). The reference condition (ref) is 
%         the same as ITE/BB. Error bars show SEM. 
%
%     `hartmann1996`    
%         Model experiments from Hartmann & Wittenberg (1996; Fig.7-8):
%         1st panel: Synthesis of zero-ILD signals. Only the harmonics 
%             from 1 to nprime had zero interaural level difference; 
%             harmonics above nprime retained the amplitudes of the baseline  
%             synthesis. Externalization scores as a function of the boundary  
%             harmonic number nprime. Fundamental frequency of 125 Hz.
%         2nd panel: Synthesis of signals to test the ISLD hypothesis. 
%             Harmonics at and below the boundary retained only the interaural 
%             spectral level differences of the baseline synthesis. Higher 
%             harmonics retained left and right baseline harmonic levels. 
%             Externalization scores as a function of the boundary
%             frequency.
%
%     `hassager2016`    
%         Model experiments from Hassager et al. (2016; Fig.6):
%         The mean of the seven listeners perceived sound source 
%         location (black) as a function of the bandwidth factor 
%         and the corresponding model predictions (colored). 
%         The model predictions have been shifted slightly to the right 
%         for a better visual interpretation. The error bars are one 
%         standard error of the mean.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2017
%
%   3) Statistics Toolbox for Matlab (for some of the figures)
%
%   Examples:
%   ---------
%
%   To display results for Fig.6 from Hassager et al. (2016) use :::
%
%     exp_baumgartner2017('hassager2016');
%
%   References: Hassager et al. (JASA 2016)
   
% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

definput.import={'amtcache'};
definput.flags.type = {'missingflag','boyd2012','hartmann1996','hassager2016'};
definput.flags.quickCheck = {'','quickCheck'};
definput.keyvals.Sintra = 2;
definput.keyvals.Sinter = 1;

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

%% Hassager et al. (2016)
if flags.do_hassager2016
  azi = [0,50];
  flp = 6000; % lowpass filtered noise stimuli
  
  Pext_A = data_hassager2016;
  B = Pext_A.B;
  
  fncache = ['hassager2016_Sintra',num2str(kv.Sintra*100,'%i'),'_Sinter',num2str(kv.Sinter*100,'%i')];
  Pext = amtcache('get',fncache,flags.cachemode);
  if isempty(Pext)
    
    data = data_baumgartner2017;
    if flags.do_quickCheck
      data = data(1:5);
    end

  %   fs = data(1).Obj.Data.SamplingRate;
  %   in = noise(fs/2,1,'white');
  %   [b,a]=butter(10,6000/fs,'low');
  %   in = filter(b,a,in);

    Pext = nan(length(B),length(data),length(azi));
    Pext = {Pext,Pext};
  %   Lat = Pext;
    for isubj = 1:length(data)
      Obj = data(isubj).Obj;
      for iazi = 1:length(azi)
  % figure;
    %     templateSound = SOFAspat(in,Obj,azi(iazi),0);
        idazi = Obj.SourcePosition(:,1) == azi(iazi) & Obj.SourcePosition(:,2) == 0;
        template = squeeze(shiftdim(Obj.Data.IR(idazi,:,:),2));
        for iB = 1:length(B)
          disp(iB)
          if isnan(B(iB))
            target = template;
          else
            Obj_tar = hassager2016spectralsmoothing(Obj,B(iB));
            target = squeeze(shiftdim(Obj_tar.Data.IR(idazi,:,:),2));
          end
  % [tarmp,fc] = baumgartner2014spectralanalysis(target);
  % subplot(1,2,1); hold on
  % semilogx(fc,tarmp(:,1))
  % % plotfftreal(iB*db2mag(10)*fftreal(target(:,1)),fs,'flog'); 
  % subplot(1,2,2); hold on
  % semilogx(fc,tarmp(:,2))
  % % plotfftreal(iB*db2mag(10)*fftreal(target(:,2)),fs,'flog'); 
          Pext{1}(iB,isubj,iazi) = baumgartner2017(target,template,'S',kv.Sinter,'flow',100,'fhigh',flp,'interaural'); % Obj instead of single template
          Pext{2}(iB,isubj,iazi) = baumgartner2017(target,template,'S',kv.Sintra,'flow',100,'fhigh',flp,'lat',azi(iazi));
        end
      end
      amtdisp([num2str(isubj),' of ',num2str(length(data)),' subjects completed.'],'progress')
    end
    
    if not(flags.do_quickCheck)
      amtcache('set',fncache,Pext);
    end
  end
  
  %% Plot
  BplotTicks = logspace(log10(0.25),log10(64),9);
  BplotTicks = round(BplotTicks*100)/100;
  BplotStr = ['Ref.';num2str(BplotTicks(:))];
  BplotTicks = [BplotTicks(1)/2,BplotTicks];
  B(1) = B(2)/2;
  Ns = size(Pext{1},2);
  
  figure 
  symb = {'-s','-d'};
  for iazi = 1:length(azi)
    subplot(1,2,iazi)
    h(3) = plot(B,Pext_A.rating(:,iazi),'-ko');
    hold on
    for m = 1:2
      h(m) = errorbar(B,mean(Pext{m}(:,:,iazi),2),std(Pext{m}(:,:,iazi),0,2)/sqrt(Ns),symb{m});
    end
    set(h,'MarkerFaceColor','w')
    set(gca,'XTick',BplotTicks,'XTickLabel',BplotStr,'XScale','log')
    axis([BplotTicks(1)/1.5,BplotTicks(end)*1.5,0.8,5.2])
    xlabel('Bandwidth Factor [ERB]')
    ylabel('Mean Externalization Rating')
    title([num2str(azi(iazi)),'\circ'])
  end
  leg = legend({'Actual','Interaural','Monaural'},'Location','southwest');
  set(leg,'Box','off')
  
end

%% Hartmann & Wittenberg (1996)
if flags.do_hartmann1996
  
  exp = {'ILD','ISLD'};
  nprime = [0,1,8,14,19,22,25,38];
  
  fncache = ['hartmann1996_Sintra',num2str(kv.Sintra*100,'%i'),'_Sinter',num2str(kv.Sinter*100,'%i')];
  Pext = amtcache('get',fncache,flags.cachemode);
  if isempty(Pext)
    azi = -37;
    data = data_baumgartner2017;
    if flags.do_quickCheck
      data = data(1:5);
    end
    Pext{1} = nan(length(data),length(nprime),length(exp));
    Pext{2} = nan(length(data),length(nprime),length(exp));
    for isub = 1:length(data)
      Obj = data(isub).Obj;
      template = hartmann1996siggen(0,'Obj',Obj,'dur',0.1);
      for ee = 1:length(exp)
        for nn = 1:length(nprime)
          target = hartmann1996siggen(nprime(nn),'Obj',Obj,'dur',0.1,exp{ee});
          Pext{1}(isub,nn,ee) = baumgartner2017(target,template,'c1',3,'c2',0,'S',kv.Sinter,'flow',100,'fhigh',6000,'interaural');
          Pext{2}(isub,nn,ee) = baumgartner2017(target,template,'c1',3,'c2',0,'S',kv.Sintra,'flow',100,'fhigh',6000,'lat',azi);
        end
      end
      amtdisp([num2str(isub),' of ',num2str(length(data)),' subjects completed.'],'progress')
    end
    if not(flags.do_quickCheck)
      amtcache('set',fncache,Pext);
    end
  end
  
  Ns = size(Pext{1},1);
  
  figure
  for ee = 1:length(exp)
    act = data_hartmann1996(exp{ee});
    subplot(1,2,ee)
    h(1) = plot(act.avg.nprime,act.avg.Escore,'-ko');
    hold on
    h(2) = errorbar(nprime,mean(Pext{1}(:,:,ee)),std(Pext{1}(:,:,ee))/sqrt(Ns),'-s');
    h(3) = errorbar(nprime,mean(Pext{2}(:,:,ee)),std(Pext{2}(:,:,ee))/sqrt(Ns),'-d');
    set(h,'MarkerFaceColor','w')
    if ee == 1
      xlabel('n\prime (Highest harmonic with ILD = 0)')
    else
      xlabel('n\prime (Highest harmonic with altered amplitudes)')
    end
    ylabel('Externalization score')
    title(exp{ee})
    axis([0,39,-0.1,3.1])
    if ee == 1
      leg = legend('Actual','Interaural cues','Monaural cues','Location','south');
      set(leg,'Box','off')
    end
  end
end

%% Boyd et al. (2012)
if flags.do_boyd2012
  flp = 6500; % Low-pass cut-off frequency
  ele = 0; 
  azi = -30;

  
  fncache = ['boyd2012_Sintra',num2str(kv.Sintra*100,'%i'),'_Sinter',num2str(kv.Sinter*100,'%i')];
  E = amtcache('get',fncache,flags.cachemode);
  if isempty(E)
    
    data = data_boyd2012;
    Eboyd = cat(3,[data.ITE.BB(:),data.BTE.BB(:)],[data.ITE.LP(:),data.BTE.LP(:)]);
    mix = data.mix/100;

    %%
    sig = noise(50e3,1,'pink');
    BTE = data_baumgartner2017('BTE');
    subjects = data_baumgartner2017;
    if flags.do_quickCheck
      subjects = subjects(1:5);
    end
    fs = subjects(1).Obj.Data.SamplingRate;

    Ebaum = nan([size(Eboyd),length(subjects)]);
    Ehass = Ebaum;
    for isub = 1:length(subjects)
      ITE = subjects(isub).Obj;

      %%
      stim{1,1} = SOFAspat(sig,ITE,azi,ele);
      stim{2,1} = SOFAspat(sig,BTE.Obj,azi,ele);
      
      % Head-absent condition
      itd = abs(round(3*0.08/343*sin(deg2rad(azi))*fs)); % |ITD| in samples
      stim{3,1} = [ [zeros(itd,1);sig(:)] , [sig(:);zeros(itd,1)] ];
      stim{3,1} = [stim{3,1};zeros(length(stim{1,1})-length(stim{3,1}),2)];
      if azi > 0
        stim{3,1} = filplr(stim{3,1});
      end
      SPL = mean(dbspl(stim{1,1}));
      stim{3,1} = setdbspl(stim{3,1},SPL);
      
%       HA = baumgartner2017flattenhrtfmag(ITE,0,50,18e3);
%       stim{3,1} = SOFAspat(sig,HA,azi,ele);
%       % remove ILD in HA condition
%       SPL = mean(dbspl(SOFAspat(sig,HA,0,0)));
%       stim{3,1} = setdbspl(stim{3,1},SPL);

      template = stim{1,1};

      %%
      % lookup = itd2anglelookuptable(ITE,ITE.Data.SamplingRate,'dietz2011');
      % [phi,phi_std,itd,ild,cfreqs] = wierstorf2013estimateazimuth(...
      %     stim{1,1},lookup,'fs',sig.fs,'dietz2011','rms_weighting');
      % % ild = mean(ild);

      %% Low pass
      [b,a]=butter(10,2*flp/fs,'low');
      for istim = 1:length(stim)
        stim{istim,2} = filter(b,a,stim{istim,1});
      end

      %%
      temSPL = dbspl(template);
      target = cell(length(mix),2,2);
      fhigh = [16e3,flp];
      for c = 1:2
        for lp = 1:2
          for m = 1:length(mix)
            target{m,c,lp} = mix(m)*stim{c,lp} + (1-mix(m))*stim{3,lp};
            for ch = 1:2
              target{m,c,lp}(:,ch) = setdbspl(target{m,c,lp}(:,ch),temSPL(ch));
            end
            Ebaum(m,c,lp,isub) =  baumgartner2017( target{m,c,lp},template,...
              'S',kv.Sintra,'flow',100,'c1',100,'c2',0,'fhigh',fhigh(1),'lat',azi);
            Ehass(m,c,lp,isub) =  baumgartner2017( target{m,c,lp},template,...
              'S',kv.Sinter,'flow',100,'c1',100,'c2',0 ,'fhigh',fhigh(1),'interaural');
          end
        end
      end
      amtdisp([num2str(isub),' of ',num2str(length(subjects)),' subjects completed.'],'progress')

    end

    seEbaum = std(Ebaum,0,4);
    Ebaum = mean(Ebaum,4);
    seEhass = std(Ehass,0,4);
    Ehass = mean(Ehass,4);
    seEboyd = 0*seEhass;

    E.m = {Eboyd,Ehass,Ebaum};
    E.se = {seEboyd,seEhass,seEbaum};
    
    if not(flags.do_quickCheck)
      amtcache('set',fncache,E);
    end
  end

  %% Plot
  figure
  dataLbl = {'Actual','Interaural','Monaural'};
  symb = {'ko','s','d'};
  condLbl = {'ITE / BB','BTE / BB','ITE / LP','BTE / LP'};
  hax = tight_subplot(1,4,0,[.15,.1],[.1,.05]);
  for cc = 1:length(condLbl)
    axes(hax(cc))
    for ee = 1:length(E.m)
      h = errorbar(data.mix,E.m{ee}(:,cc),E.se{ee}(:,cc),['-',symb{ee}]);
      set(h,'MarkerFaceColor','w')
      hold on
    end
    set(gca,'XDir','reverse')
    xlabel('Mix (%)')
    if cc == 1
      ylabel('Externalization score')
    else
      set(gca,'YTickLabel',[])
    end
    title(condLbl{cc})
    axis([-20,120,-5,105])
    if cc == 1
      leg = legend(dataLbl);
      set(leg,'Box','off','Location','south')
    end
  end
%   set(leg,'Location','eastoutside','Position',get(leg,'Position')+[.1,.2,0,0])
end







%%%%%%%%%% INTERNAL FUNCTIONS FOR VISUALIZATION %%%%%%%%%%
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

end