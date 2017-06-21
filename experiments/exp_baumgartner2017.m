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

definput.flags.type = {'missingflag','boyd2012','hartmann1996','hassager2016'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

%% Hassager et al. (2016)
if flags.do_hassager2016
  azi = [0,50];
  
  Pext_A = data_hassager2016;
  B = Pext_A.B;
  
  data = data_baumgartner2017;
%   data = data(1:5);
  fs = data(1).Obj.Data.SamplingRate;
%   Obj = SOFAload('BRIR_AllAbsorbers_OffCentre_Emitters1to64_front.sofa');
%   Obj = SOFAload('dtf b_nh10.sofa');
  in = noise(fs/2,1,'white');
  [b,a]=butter(10,6000/fs,'low');
  in = filter(b,a,in);
  
  Pext = nan(length(B),length(data),length(azi));
  Pext_B = Pext;
  Lat = Pext;
  for isubj = 1:length(data)
    Obj = data(isubj).Obj;
    for iazi = 1:length(azi)
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
  % plotfftreal(iB*db2mag(10)*fftreal(target(:,1)),fs,'flog'); hold on
%         Pext(iB,isubj,iazi) = hassager2016(target,template,'fs',fs);
        [Pext(iB,isubj,iazi),Lat(iB,isubj,iazi)] = baumgartner2017(target,template,'flow',50,'fhigh',6000,'lat',azi(iazi),'interaural'); % Obj instead of single template
        Pext_B(iB,isubj,iazi) = baumgartner2017(target,template,'S',5,'flow',100,'fhigh',18e3,'lat',azi(iazi));
      end
%       Pext_B(:,isubj,iazi) = Pext_B(:,isubj,iazi)*4+1; % scale to be within 1 to 5   %/Pext_B(1,isubj,iazi)
    end
    disp([num2str(isubj),' of ',num2str(length(data)),' subjects completed.'])
  end
  
  BplotTicks = logspace(log10(0.25),log10(64),9);
  BplotTicks = round(BplotTicks*100)/100;
  BplotStr = ['Ref.';num2str(BplotTicks(:))];
  BplotTicks = [BplotTicks(1)/2,BplotTicks];

  %%
  B(1) = B(2)/2; % just for
  figure 
  for iazi = 1:length(azi)
    subplot(1,2,iazi)
    plot(B,Pext_A.rating(:,iazi),'k')
    hold on
    h(1) = errorbar(B,mean(Pext(:,:,iazi),2),std(Pext(:,:,iazi),0,2)/sqrt(length(data)),'-s');
    h(2) = errorbar(B,mean(Pext_B(:,:,iazi),2),std(Pext_B(:,:,iazi),0,2)/sqrt(length(data)),'-d');
    set(h,'MarkerFaceColor','w')
    set(gca,'XTick',BplotTicks,'XTickLabel',BplotStr,'XScale','log')
    axis([BplotTicks(1)/1.5,BplotTicks(end)*1.5,0.8,5.2])
    xlabel('Bandwidth Factor [ERB]')
    ylabel('Mean Externalization Rating')
    title([num2str(azi(iazi)),'\circ'])
  end
  leg = legend({'Actual','Interaural','Monaural'},'Location','southwest');
  set(leg,'Box','off')
%   RB_print(gcf,[14,6],'exp_hassager2016')
  
end

%% Hartmann & Wittenberg (1996)
if flags.do_hartmann1996
  
  data = data_baumgartner2017;
%   data = data(1:5);
  nprime = [0,1,8,14,19,22,25,38];
  exp = {'ILD','ISLD'};
  Pext{1} = nan(length(data),length(nprime),length(exp));
  Pext{2} = nan(length(data),length(nprime),length(exp));
  for isub = 1:length(data)
    Obj = data(isub).Obj;
    template = sig_hartmann1996(0,'Obj',Obj,'dur',0.1);
    for ee = 1:length(exp)
      for ii = 1:length(nprime)
        target = sig_hartmann1996(nprime(ii),'Obj',Obj,'dur',0.1,exp{ee});
%         Pext{1}(isub,ii,ee) = hassager2016(target,template,'c1',3,'c2',0,'flow',100);
  %       E(ii,ee) = 3*baumgartner2015externalization(target,template,'S',0.4);
        Pext{1}(isub,ii,ee) = baumgartner2017(target,template,'c1',3,'c2',0,'S',1,'flow',100,'fhigh',6000,'interaural');
        Pext{2}(isub,ii,ee) = baumgartner2017(target,template,'c1',3,'c2',0,'S',5,'flow',100,'fhigh',18e3);
      end
%       Pext{2}(isub,:,ee) = 3*Pext{2}(isub,:,ee)/Pext{2}(isub,1,ee);
    end
    disp([num2str(isub),' of ',num2str(length(data)),' subjects completed.'])
  end
  
  figure
  for ee = 1:length(exp)
    act = data_hartmann1996(exp{ee});
    subplot(1,2,ee)
    plot(act.avg.nprime,act.avg.Escore,'k');
    hold on
    h(1) = errorbar(nprime,mean(Pext{1}(:,:,ee)),std(Pext{1}(:,:,ee))/sqrt(length(data)),'-s');
    h(2) = errorbar(nprime,mean(Pext{2}(:,:,ee)),std(Pext{2}(:,:,ee))/sqrt(length(data)),'-d');
    set(h,'MarkerFaceColor','w')
    if ee == 1
      xlabel('n\prime (Highest harmonic with ILD = 0)')
    else
      xlabel('n\prime (Highest harmonic with altered amplitudes)')
    end
    ylabel('Externalization score')
    title(exp{ee})
    axis([0,39,-0.1,3.1])
  end
  leg = legend('Actual','Interaural cues','Monaural cues','Location','east');
  set(leg,'Box','off')
%   RB_print(gcf,[12,6],'exp_hassager2016_hartmann1996')
end

%% Boyd et al. (2012)
if flags.do_boyd2012
  flp = 6500; % Low-pass cut-off frequency
  ele = 0; 
  azi = -30;

  data = data_boyd2012;
  Eboyd = cat(3,[data.ITE.BB(:),data.BTE.BB(:)],[data.ITE.LP(:),data.BTE.LP(:)]);
  mix = data.mix/100;

  %%
  sig.scare = noise(50e3,1,'pink');
  BTE = data_baumgartner2017('BTE');
  subjects = data_baumgartner2017;
  % subjects = subjects(1:3);
  fs = subjects(1).Obj.Data.SamplingRate;

  Ebaum = nan([size(Eboyd),length(subjects)]);
  Ehass = Ebaum;
  for isub = 1:length(subjects)
    ITE = subjects(isub).Obj;
    HA = baumgartner2017flattenhrtfmag(ITE,0,50,18e3);

    %%
    stim{1,1} = SOFAspat(sig.scare,ITE,azi,ele);
    stim{2,1} = SOFAspat(sig.scare,BTE,azi,ele);
    stim{3,1} = SOFAspat(sig.scare,HA,azi,ele);
    % remove ILD in HA condition
    SPL = mean(dbspl(SOFAspat(sig.scare,HA,0,0)));
    stim{3,1} = setdbspl(stim{3,1},SPL);

    template = stim{1,1};%SOFAspat(sig.noise,ITE,azi,ele);%stim{1,1};
    % template = permute(stim{1,1},[1,3,2]);
    % template = SOFAspat(noise(sig.fs/10,1,'white'),ITE,azi,ele);

    %%
    % lookup = itd2anglelookuptable(ITE,ITE.Data.SamplingRate,'dietz2011');
    % [phi,phi_std,itd,ild,cfreqs] = wierstorf2013estimateazimuth(...
    %     stim{1,1},lookup,'fs',sig.fs,'dietz2011','rms_weighting');
    % % ild = mean(ild);

    %% Low pass
    [b,a]=butter(10,2*flp/fs,'low');
    for ii = 1:length(stim)
      stim{ii,2} = filter(b,a,stim{ii,1});
    end

    %%

    target = cell(length(mix),2,2);
    fhigh = [16e3,flp];
    for c = 1:2
      for lp = 1:2
        for m = 1:length(mix)
          target{m,c,lp} = mix(m)*stim{c,lp} + (1-mix(m))*stim{3,lp};
          Ebaum(m,c,lp,isub) =  baumgartner2017( target{m,c,lp},template,...
            'S',5,'flow',100,'c1',100,'c2',0,'fhigh',fhigh(1));
          Ehass(m,c,lp,isub) =  baumgartner2017( target{m,c,lp},template,...
            'S',1,'flow',100,'c1',100,'c2',0 ,'fhigh',fhigh(1),'interaural');
    %       Ehass(m,c,lp,isub) =  hassager2016( target{m,c,lp},template,...
    %                                     'fhigh',fhigh(lp),'c1',100,'c2',0 );
        end
      end
    end
    amtdisp([num2str(isub),' of ',num2str(length(subjects)),' subjects completed.'])

  end

  seEbaum = std(Ebaum,0,4);
  Ebaum = mean(Ebaum,4);
  seEhass = std(Ehass,0,4);
  Ehass = mean(Ehass,4);
  seEboyd = 0*seEhass;

  %%
  figure
  E = {Eboyd,Ehass,Ebaum};
  seE = {seEboyd,seEhass,seEbaum};
  dataLbl = {'Actual','Interaural','Monaural'};
  symb = {'ko','s','d'};
  condLbl = {'ITE / BB','BTE / BB','ITE / LP','BTE / LP'};
  hax = tight_subplot(1,4,0,[.15,.1],[.1,.15]);
  for ii = 1:length(condLbl)
    axes(hax(ii))
    for ee = 1:length(E)
      h = errorbar(data.mix,E{ee}(:,ii),seE{ee}(:,ii),['-',symb{ee}]);
      set(h,'MarkerFaceColor','w')
      hold on
    end
    set(gca,'XDir','reverse')
    xlabel('Mix (%)')
    if ii == 1
      ylabel('Externalization score')
    else
      set(gca,'YTickLabel',[])
    end
    title(condLbl{ii})
    axis([-20,120,-5,105])
  end
  leg = legend(dataLbl,'Location','eastoutside');
  set(leg,'Position',get(leg,'Position')+[.1,.2,0,0])

  % RB_print(gcf,[13,6],'exp_baumgartner2017_boyd2012')
end

end