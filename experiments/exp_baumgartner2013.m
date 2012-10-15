function varargout=exp_baumgartner2013(varargin)
%EXP_BAUMGARTNER2013   Figures from Baumgartner et al. (2013)
%   Usage: output = exp_baumgartner2013(flag)
%
%   `exp_baumgartner2013(flags,... )` reproduces figures of the book 
%   chapter from Baumgartner et al. (2013).
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'fig12'   Reproduce Fig. 12:
%               DTFs of median SP (left ear). 
%               Left: NH12. Center: NH58. Right: NH33. 
%               Brightness: Spectral magnitude.
%
%     'fig13'   Reproduce Fig. 13:
%               Listener's localization performance for non-individualized 
%               versus individualized DTFs. 
%               Bars: Individualized DTFs. 
%               Circles: Non-individualized DTFs averaged over 16 DTF sets. 
%               Error bars: ±1 standard deviation of the average. 
%               Horizontal line: Chance performance corresponding to 
%               guessing the direction of the sound.
%
%     'fig14'   Reproduce Fig. 14:
%               Bars: Increase in localization errors when listening to a 
%               binaural recording created with the DTFs from NH58. 
%               Asterisk: Difference between actual (experiment) and 
%               predicted performance when listening to individualized DTFs. 
%               Horizontal lines: chance performance corresponding to 
%               guessing the direction of the sound.
%
%     'fig15'   Reproduce Fig. 15:
%               Bars: Localization performance of the pool listening to 
%               selected DTFs. 
%               Circles: DTFs from NH12. 
%               Squares: DTFs from NH58. 
%               Triangles: DTFs from NH33.
%
%     'fig18'   Reproduce Fig. 18:
%               Predicted response probabilities (PMVs) as a function of 
%               the amplitude panning ratio between a rear and a front 
%               loudspeaker in the same SP (lateral angle: 30°). 
%               Left: Results for NH12. 
%               Center: Results for NH15. 
%               Right: Results for our pool of listeners.
%               Circle: Maximum of a PMV. 
%               Panning ratio of 0: Only front loudspeaker active. 
%               Panning ratio of 1: Only rear loudspeaker active. 
%
%   See also: baumgartner2013, data_baumgartner2013
%
%   Examples:
%   ---------
%
%   To display Fig. 13 use :::
%
%     exp_baumgartner2013('fig13');
%
%   References:baumgartner2013,lopezpoveda2001hnc pralong1996role goode1994nkf
  
% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

%% ------ Check input options --------------------------------------------

  definput.flags.type = {'missingflag',...
    'fig12','fig13','fig14','fig15',... % Ch. 4.1 (binaural rec.)
    'fig18','fig22','fig23'}; % Ch. 4
  definput.flags.plot = {'plot','noplot'};

  % Parse input options
  [flags,kv]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% Figure 12

if flags.do_fig12
  
  s = data_baumgartner2013('pool');
  
  s = s([2 1 14]);  % NH12, NH58, NH33
  
  nfft = 2^12;
  range = 50;   % dB
  ilat = 2;     % median sagittal plane
  fhigh = 18e3; % Hz
  da = 3;       % white frame for ticks (angles)
  df = 0.1;     % white frame for ticks (freq)
  
  figure
  for ll = 1:length(s)
  
    f = 0:s(1).fs/nfft:s(ll).fs/2;
    f = f/1e3; % in kHz

    mag = db(abs(fft(s(ll).spdtfs{ilat},nfft)));
    mag = mag(1:nfft/2+1,:,:);
    maxmag = max(max(max(mag)));
    mag = mag - maxmag;       % normalize to 0dB
    maxmag = 0;
    mag(mag < maxmag-range) = maxmag-range;  % limit dynamic range
    crange = [floor(maxmag-range), ceil(maxmag)];% range of color code

    idf = f<fhigh/1e3;
    magl = mag(idf,:,1);
    fp = f(idf);
    polp = s(ll).polangs{ilat};

    subplot(1,3,ll)
    h = pcolor(fp,polp,magl');
    caxis(crange)
    shading interp
    colormap('bone')
    xlabel('Frequency (kHz)')
    ylabel('Polar Angle (°)')
    set(gca,...
        'XLim',[fp(1)-df,fp(end)+df],...
        'YLim',[polp(1)-da,polp(end)+da],...    
        'XTick',1:2:17,...
        'YTick',-30:30:polp(end),...
        'XMinorTick','on',...
        'YMinorTick','on'...
        )
    
  end
  
  if nargout == 1
    varargout{1} = s;
  end
  
end


%% Figures 13-15 (non-individual binaural recordings)

if flags.do_fig13 || flags.do_fig14 || flags.do_fig15
  
  latdivision = [-20,0,20];  % lateral center angles of SPs
  
  s = data_baumgartner2013('pool');
  ns = length(s);

  % DTFs of the SPs
  for ll = 1:ns
    for ii = 1:length(latdivision)
      
      s(ll).latang{ii} = latdivision(ii);
      s(ll).polangs{ii} = [];
      s(ll).spdtfs{ii} = [];
      [s(ll).spdtfs{ii},s(ll).polangs{ii}] = extractsp(...
        s(ll).latang{ii},s(ll).dtfs,s(ll).pos);
      
    end
  end
  
  h = waitbar(0,'Please wait a little!');
  qe = zeros(ns,ns,length(latdivision)); % init QEs
  pe = qe;                               % init PEs
  for ll = 1:ns    % listener
      for jj = 1:ns    % ears
          for ii = 1:length(latdivision) % SPs
            
            s(ll).p{jj,ii} = [];
            s(ll).respangs{ii} = [];
            [s(ll).p{jj,ii},s(ll).respangs{ii}] = baumgartner2013(...
                s(jj).spdtfs{ii},s(ll).spdtfs{ii},s(ll).fs,...
                'u',s(ll).u,'lat',s(ll).latang{ii},...
                'polsamp',s(ll).polangs{ii});

            [ qe(ll,jj,ii),pe(ll,jj,ii) ] = pmv2ppp( ...
                s(ll).p{jj,ii} , s(jj).polangs{ii} , s(ll).respangs{ii});

          end
      end
      waitbar(ll/ns,h)

  end
  close(h)

  qe = mean(qe,3);
  pe = mean(pe,3);
  for ll = 1:length(s)
    s(ll).pe = pe(ll,:);
    s(ll).qe = qe(ll,:);
  end 
  
  if nargout == 1
    varargout{1} = s; % SAME for Fig. 13-15!!!!
  end
  
  if flags.do_fig13 % -----------------------------------------------------
  
    if flags.do_plot

      markersize = 4;

      [qechance,pechance] = pmv2ppp(ones(49,44));
      
      % IDs for XTick of plots
      for ll = 1:ns
          nhs(ll,:) = s(ll).id;
      end

      % Polar error

      figure
      subplot(121)
      peinc = zeros(ns,ns-1);
      for l = 1:ns
          peinc(l,:) = pe(l,[1:l-1,l+1:ns]);
      end

      mpeinc = mean(peinc,2);
      speinc = std(peinc,0,2);

      h = bar(diag(pe));
      set(h,'FaceColor','white')

      hold on

      errorbar(mpeinc,speinc,'ko',...
        'MarkerSize',markersize,'MarkerFaceColor','w');

      plot(0:ns+1,pechance*ones(ns+2,1),'k--')

      ylabel('Local Polar RMS Error (°)')
      xlabel({'NH'})

      set(gca,'YLim',[20.1 48],'XLim',[0 ns+1],...
        'XTick',1:ns,'XTickLabel',nhs(:,3:4),...
          'YMinorTick','on')

      % Quadrant error

      subplot(122)
      qeinc = zeros(ns,ns-1);
      for l = 1:ns
          qeinc(l,:) = qe(l,[1:l-1,l+1:ns]);
      end

      mqeinc = mean(qeinc,2);
      sqeinc = std(qeinc,0,2);

      h = bar(diag(qe));
      set(h,'FaceColor','white')

      hold on
      errorbar(mqeinc,sqeinc,'ko',...
        'MarkerSize',markersize,'MarkerFaceColor','w');

      plot(0:ns+1,qechance*ones(ns+2,1),'k--')

      set(gca,'YLim',[floor(min(diag(qe))-1) ceil(qechance+1)],...
        'XLim',[0 ns+1],'XTick',1:ns,'XTickLabel',nhs(:,3:4),...
        'YMinorTick','on')

      ylabel('Quadrant Error (%)')
      xlabel({'NH'})

    end
    
  elseif flags.do_fig14 % -------------------------------------------------
    
    ears = 1;  % NH58
    
    if flags.do_plot
      
      [qechance,pechance] = pmv2ppp(ones(49,44));
      
      % IDs for XTick of plots
      for ll = 1:ns
          nhs(ll,:) = s(ll).id;
      end
      
      figure
      
      % Polar error
      subplot(121)
      h = bar(pe(:,ears)-diag(pe));
      hold on
%       plot(ears,exp.s(ears).pe_exp-pe(ears,ears),'k*')

      dx = 0.42;
      peincchance = pechance-diag(pe);
      for ll = 1:length(s)
          plot(ll+[-dx,+dx],ones(2,1)*peincchance(ll),'k--')
      end

      ylabel('Increase in Local Polar RMS Error (°)')
      xlabel('NH')
      set(h,'FaceColor','white')
      set(gca,'XLim',[0 length(s)+1],...
        'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
        'YLim',[-1 21.9],'YMinorTick','on')

      % Quadrant error
      subplot(122)
      h = bar(qe(:,ears)-diag(qe));
      hold on
%       plot(ears,exp.s(ears).qe_exp-qe(ears,ears),'k*')

      dx = 0.4;
      qeincchance = qechance-diag(qe);
      for ll = 1:length(s)
          plot(ll+[-dx,+dx],ones(2,1)*qeincchance(ll),'k--')
      end

      set(h,'FaceColor','white')
      set(gca,'XLim',[0 length(s)+1],...
        'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
        'YLim',[-5 18],'YMinorTick','on')
      xlabel('NH')
      ylabel('Increase in Quadrant Error (%)')
      
    end
    
  elseif flags.do_fig15
    
    pepoor = pe(:,14); % NH33 (14)
    pemed = pe(:,1);   % NH58 (1)
    pegood = pe(:,2);  % NH12 (2)

    qepoor = qe(:,14);
    qemed = qe(:,1);
    qegood = qe(:,2);
    
    if flags.do_plot
      
      markersize = 4;

      [qechance,pechance] = pmv2ppp(ones(49,44));
      
      % IDs for XTick of plots
      for ll = 1:ns
          nhs(ll,:) = s(ll).id;
      end
      
      figure
      
      % Polar error
      subplot(121)
      h = bar(diag(pe));
      hold on
      plot(pegood,'ko','MarkerSize',markersize,'MarkerFaceColor','k')
      plot(pemed, 'ks','MarkerSize',markersize,'MarkerFaceColor',0.5*ones(3,1))
      plot(pepoor,'kv','MarkerSize',markersize,'MarkerFaceColor','w')
      plot(0:length(s)+1,pechance*ones(length(s)+2,1),'k--')

      set(gca,'XLim',[0 length(s)+1],'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
              'YLim',[20.1 48],'YMinorTick','on')
      set(h,'FaceColor','white')%,'BarWidth',0.5)

      ylabel('Local Polar RMS Error (°)')
      xlabel('NH')

      % Quadrant error
      subplot(122)
      h = bar(diag(qe));
      hold on
      plot(qegood,'ko','MarkerSize',markersize,'MarkerFaceColor','k')
      plot(qemed,'ks','MarkerSize',markersize,'MarkerFaceColor',0.5*ones(3,1))
      plot(qepoor,'kv','MarkerSize',markersize,'MarkerFaceColor','w')
      plot(0:length(s)+1,qechance*ones(length(s)+2,1),'k--')

      set(gca,'XLim',[0 length(s)+1],'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
          'YLim',[1 44],'YMinorTick','on')
      set(h,'FaceColor','white')%,'BarWidth',0.5)

      ylabel('Quadrant Error (%)')
      xlabel('NH')
      
    end
    
  end
    
end


%% Figures 18, 22, 23 (Amplitude panning between loudspeakers)

if flags.do_fig18 || flags.do_fig22 || flags.do_fig23

  % Coordinates (azi,ele) of loudspeakers
  if flags.do_fig18
    LSPsetup{1} = [ 30,0 ; 150,0 ];                  % lat: 30°
  elseif flags.do_fig22
    LSPsetup{1} = [ 30,0 ; 110,0 ];                  % 5.1
    LSPsetup{2} = [ 30,0 ; 45,45 ; 110,0 ];          % 10.2
    LSPsetup{3} = [ 30,0 ; 30,30 ; 110,30 ; 110,0 ]; % 9.1
  elseif flags.do_fig23
    LSPsetup{1} = [ 30,0 ; 30,30 ; 110,30 ; 110,0 ]; % 9.1
    LSPsetup{2} = [ 30,0 ; 30,30 ; 135,30 ; 110,0 ]; % 9.1 (prop. LSH)
    LSPsetup{3} = [ 30,0 ; 30,30 ; 135,30 ; 135,0 ]; % 9.1 (prop. LSH&LS)
  end
  
  dpol = 5;    % step size in deg

  s = data_baumgartner2013('pool');
  ns = length(s);
  
  for ss = 1:length(LSPsetup)
  
    LSPsph = LSPsetup{ss};
    nLSP = size(LSPsph,1);
    LSPhp = zeros(nLSP,2);
    [LSPhp(:,1),LSPhp(:,2)] = sph2horpolar(LSPsph(:,1),LSPsph(:,2));
    LSPcart = zeros(nLSP,3);
    [LSPcart(:,1),LSPcart(:,2),LSPcart(:,3)] = ...
          sph2cart(LSPsph(:,1),LSPsph(:,2),ones(nLSP,1));

    D = LSPhp(1,2):dpol:LSPhp(nLSP,2);     % Desired polar angle    

    w = waitbar(0,'Please wait a little!');
    for ll = 1:ns

      [poscart(:,1),poscart(:,2),poscart(:,3)] = ...
          sph2cart(s(ll).pos(:,1),s(ll).pos(:,2),ones(size(s(ll).pos,1),1));

      % DTF indices of LSPs
      idLSP = zeros(nLSP,1);
      for ii = 1:nLSP
        [~,idLSP(ii)] = min( dist(poscart,LSPcart(ii,:)') );
      end
      clear poscart

      Lp = 1;  % LSP pair index
      for d = 1:length(D)

        if D(d) > LSPhp(Lp+1,2)  % Switch to next LSP pair?
          Lp = Lp + 1;
        end
        pol1 = LSPhp(Lp,2);
        pol2 = LSPhp(Lp+1,2);

        if abs(pol1-pol2) < 180

          pol0 = mean([pol1 pol2]);   % polar angle of center
          phi0 = abs(pol1-pol2)/2;    % basis angle
          M = (1-tan(deg2rad(pol0-D(d)))/tan(deg2rad(phi0)))/2;

        else

          M = (D(d)-pol1) / (pol2-pol1);

        end

        lat1 = LSPhp(Lp,1);
        lat2 = LSPhp(Lp+1,1);
        lat0 = mean([lat1 lat2]);   % lateral angle of center
        phi0 = abs(lat1-lat2)/2;    % basis angle
        phiM = rad2deg( atan( tan(deg2rad(phi0)) * (1-2*M) ) );
        latM = lat0 - phiM;

        dtfLSP = (1-M)*s(ll).dtfs(:,idLSP(Lp),:) + M*s(ll).dtfs(:,idLSP(Lp+1),:);
        [spdtfs,polangs] = extractsp(latM,s(ll).dtfs,s(ll).pos);

        [pmv,respangs] = baumgartner2013(dtfLSP,spdtfs,s(ll).fs,...
            'lat',latM,'u',s(ll).u,...
            'polsamp',polangs);
        idP = respangs >= -30 & respangs <= 210;
        s(ll).p{ss}(:,d) = pmv(idP)/sum(pmv(idP));
        s(ll).respangs = respangs(idP);

      end

      waitbar(ll/ns,w)

    end
    close(w)

    % Pool
    ppool = zeros(size(s(1).p{ss}));
    for ll = 1:ns
        ppool = ppool + s(ll).p{ss};
    end
    s(ns+1).p{ss} = ppool/length(s);
    s(ns+1).respangs = s(1).respangs;
    s(ns+1).id = 'Pool';
  
  end
  
  % Output
  if nargout == 1
    varargout{1} = s;
  end
  
  % Plots
  if flags.do_plot
    
    if flags.do_fig18
      subject = [2 10 18];
      setup = [1 1 1];
    elseif flags.do_fig22
      subject = [18 18 18];
      setup = 1:3;
    elseif flags.do_fig23
      subject = [18 18 18];
      setup = 1:3;
    end
    
    figure
    
    for ii = 1:3
      ll = subject(ii);
      ss = setup(ii);
      subplot(1,3,ii)
    
      h = pcolor(D,s(ll).respangs,s(ll).p{ss});
      axis equal
      if flags.do_fig23 || (flags.do_fig22 && ss > 1)
        set(gca,'XLim',[D(1)-3 D(end)+3],'XMinorTick','on',...
          'XTick',0:30:180,'XTickLabel',{0:30:150,180})
        xlabel('Desired Polar Angle (°)')
      else
        set(gca,'XLim',[D(1)-3 D(end)+3],'XMinorTick','on',...
          'XTick',18:36:180,'XTickLabel',{0.1:0.2:0.9})
        xlabel('Panning Ratio')
      end
      set(gca,'YLim',[s(ll).respangs(1)-3 s(ll).respangs(end)+3],...
      	'YMinorTick','on','YTick',(0:30:180)+2.5,'YTickLabel',0:30:180)
      colormap bone
      shading flat
      caxis([0 0.1])
      ylabel('Response Angle (°)')

      [tmp,Imax] = max(s(ll).p{ss});
      hold on
      h1 = plot( D,s(ll).respangs(Imax)+2.5, 'wo');  % shadow
      set(h1,'MarkerSize',4,'MarkerFaceColor','none') 
      h2 = plot( D,s(ll).respangs(Imax)+2.5, 'ko'); 
      set(h2,'MarkerSize',3,'MarkerFaceColor','none') 
      
    end
    
  end
  
end

end