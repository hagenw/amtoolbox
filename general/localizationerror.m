function [varargout] = localizationerror(m,errorflag)
%LOCALIZATIONERROR Compute psychoacoustic performance parameters for sound
%localization experiments
%
%   Usage: [accL, precL, accP, precP, querr] = localizationerror(m)
%          [res, meta, par] = localizationerror(m,errorflag)
%
%   `localizationerror(m,errorflag)` returns psychoacoustic performance 
%   parameters for experimental response patterns .
%   `[res, meta, par] = localizationerror(m,errorflag)` calculates error 
%   metric *res* given by the string *errorflag*. *meta* contains additional 
%   metadata, *par* contains additional calculation results.
%
%   The errorflag may be one of:
%
%   accL:       accuracy (average) in the lateral dimension
%   precL:      precision (std. dev.) in the lateral dimension
%   accP:       accuracy (average) in the polar dimension
%   precP:      precision (circular std. dev.) in the polar dimension 
%   querr:      percentage of weighted polar errors > 45°: 
%               [querr number_of_corrects total_number] = CalcAccPrec...
%   accE:       accuracy in the elevation
%   absaccE:    absolute accuracy in the elevation
%   absaccL:    absolute lateral accuracy 
%   accabsL:    mean absolute lateral error
%   accPnoquerr:accP with quadrant errors removed
%   precE:      precision in the elevation
%   querr90:    precentage of weighted polar errors > 90°
%   precPmedian:unweighted polar precision error considering targets around
%                   the median plane only (+/-30° lateral) 
%   precPmedianlocal: basing on precPmedian, quadrant errors (polar error > 90°)
%                   are excluded. Similar to RMS local polar error (Middlebrooks, 1999)
%   precPnoquerr:precP with quadrant errors removed
%   rmsL:       lateral RMS error according to Middlebrooks (1999)
%   rmsPmedianlocal: unweighted polar RMS considering targets around
%                   the median plane only (+/-30° lateral) and excluding quadrant
%                   errors (polar error > 90°). Identical to RMS local polar
%                   error in Middlebrooks (1999).
%   rmsPmedian: unweighted polar RMS considering targets around the median plane 
%                   only (+/-30° lateral). Includes quadrant errors.
%   querrMiddlebrooks: quadrant errors as in Middlebrooks (1999). Use
%                   central positions within +/-30°, querr is when polar error is >90°
%                   [querr number_of_corrects total_number] = CalcAccPrec...
%   corrcoefL:  lateral correlation coeffficient [cc, p] = CalcAccPrec...
%   corrcoefP:  polar correlation coeffficient, Only data within lateral +/-30° 
%                   are considered. [corr_coeff, p_of_significant_corr] = CalcAccPrec...
%   SCC:        spherical/spatial correlation coefficient 
%                   (Carlile et al., 1997; Fisher et al., 1987)

%   References: FIXME

% AUTHOR : Piotr Majdak, Robert Baumgartner


if size(m,2)==1,   m=m'; end

    
  if ~exist('errorflag','var')
    nargout=6;
    accL=mean(m(:,7)-m(:,5));
    varargout{1}=accL;
    precL=sqrt(sum((m(:,7)-m(:,5)-accL).^2)/size(m,1));  
    varargout{2}=precL;
    w=0.5*cos(deg2rad(m(:,7))*2)+0.5;  
    varargout{3}=npi2pi(circmean(m(:,8)-m(:,6),w));
    [x,mag]=circmean(m(:,8)-m(:,6));
    precP=rad2deg(acos(mag)*2)/2;        
    varargout{4}=precP;
    idx=find((abs(npi2pi(m(:,8)-m(:,6))).*w)<=thr);  
    querr=100-size(idx,1)/size(m,1)*100;
    varargout{5}=querr;
  else
    if iscell(errorflag), errorflag=errorflag{1}; end
    if ~ischar(errorflag)
      error('ERR must be a string');
    end
    if size(errorflag,1)~=1, errorflag=errorflag'; end
    nargout=3;
    metadata=[];
    par=[];
    switch errorflag
      case 'accL'
        accL=mean(m(:,7)-m(:,5));  
        varargout{1}=accL;
        meta.ylabel='Lateral bias (deg)';
      case 'absaccL'
        accL=mean(m(:,7)-m(:,5));  
        varargout{1}=abs(accL);
        meta.ylabel='Lateral absolute bias (deg)';
      case 'accabsL'
        accL=mean(abs(m(:,7)-m(:,5)));  
        varargout{1}=accL;
        meta.ylabel='Lateral absolute error (deg)';
      case 'corrcoefL'
        [r,p] = corrcoef([m(:,5) m(:,7)]);
        varargout{1}=r(1,2);
        par.par1=p(1,2);
        meta.ylabel='Lateral correlation coefficient';
        meta.par1 = 'p for testing the hypothesis of no correlation';
      case 'corrcoefP'
        idx=find(abs(m(:,7))<=30);
        [r,p] = corrcoef([m(idx,6) m(idx,8)]);
        varargout{1}=r(1,2);
        par.par1=p(1,2);
        meta.ylabel='Polar correlation coefficient';
        meta.par1 = 'p for testing the hypothesis of no correlation';
      case 'precL'
        accL=mean(m(:,7)-m(:,5));  
        precL=sqrt(sum((m(:,7)-m(:,5)-accL).^2)/size(m,1));  
        varargout{1}=precL;
        meta.ylabel='Lateral precision error (deg)';
      case 'rmsL'
        idx=find(abs(m(:,7))<=60);
        m=m(idx,:);        
        rmsL=sqrt(sum((npi2pi(m(:,7)-m(:,5))).^2)/size(m,1));
        varargout{1}=rmsL;
        meta.ylabel='Lateral RMS error (deg)';
      case 'accP'
        varargout{1}=npi2pi(circmean(m(:,8)-m(:,6)));
        meta.ylabel='Polar bias (deg)';
      case 'accPnoquerr'
        w=0.5*cos(deg2rad(m(:,7))*2)+0.5;
        idx=find(w>0.5);
        m=m(idx,:); w=w(idx,:);
        idx=find((abs(npi2pi(m(:,8)-m(:,6))).*w)<=45);                  
        varargout{1}=npi2pi(circmean(m(idx,8)-m(idx,6)));     
        par.par1=length(idx);
        meta.ylabel='Local polar bias (deg)';
        meta.par1='Number of considered responses';
      case 'precP'
        w=0.5*cos(deg2rad(m(:,7))*2)+0.5;
        precP=circstd(m(:,8)-m(:,6),w);
        varargout{1}=precP;
        meta.ylabel='Polar precision (deg)';
      case 'precPnoquerr'
%         w=0.5*cos(deg2rad(m(:,7))*2)+0.5;
%         idx=find(w>0.5);
%         thr=45; % threshold for the choice: quadrant error or not        
%         m=m(idx,:); w=w(idx,:);
        range=360; % polar range of targets
        thr=45; % threshold for the choice: quadrant error or not
        latmax=0.5*180*(acos(4*thr/range-1))/pi;
        idx=find(abs(m(:,7))<latmax); % targets outside not considered
        m=m(idx,:);
        w=0.5*cos((m(:,7))*2*pi/180)+0.5;
        idx=find((abs(npi2pi(m(:,8)-m(:,6))).*w)<=thr);                  
        m=m(idx,:);
        w=0.5*cos((m(:,7))*pi/180*2)+0.5;
        precP=circstd(m(:,8)-m(:,6),w);
        varargout{1}=real(precP);        
        meta.ylabel='Local polar precision error (deg)';
      case 'querr'
        range=360; % polar range of targets
        thr=45; % threshold for the choice: quadrant error or not
        latmax=0.5*180*(acos(4*thr/range-1))/pi;
        idx=find(abs(m(:,7))<latmax); % targets outside not considered
        m=m(idx,:);
        w=0.5*cos(pi*(m(:,7))/180*2)+0.5; % weigthing of the errors
        corrects=size(find((abs(npi2pi(m(:,8)-m(:,6))).*w)<=thr),1);
        totals=size(m,1);
        chancerate=2*thr/range./w;
        wrongrate=100-mean(chancerate)*100;
        querr=100-corrects/totals*100;
        varargout{1}=querr;
        meta.ylabel='Weighted quadrant errors (%)';
        par.par1=corrects;
        meta.par1='Number of correct responses';
        par.par2=totals;
        meta.par2='Number of responses within the lateral range'; 
        par.par3=wrongrate; % valid only for targets between 0 and 180°
        meta.par3='Wrong rate, valid only for targets between 0 and 180°';
      case 'querrMiddlebrooks' % as in Middlebrooks (199); comparable to COC from Best et al. (2005)
        idx=find(abs(m(:,7))<=30);
        m=m(idx,:);
        idx=find((abs(npi2pi(m(:,8)-m(:,6))))<=90);  
        querr=100-size(idx,1)/size(m,1)*100;
        varargout{1}=querr;
        meta.ylabel='Quadrant errors (%)';
        par.par1=size(idx,1);
        meta.par1='Number of confusions';
        par.par2=size(m,1);
        meta.par2='Number of responses within the lateral range';
      case 'querrMiddlebrooksRAU' % as in Middlebrooks (1999) but calculated in rau
        idx=find(abs(m(:,7))<=30);
        m=m(idx,:);
        idx=find((abs(npi2pi(m(:,8)-m(:,6))))<=90);  
        querr=100-size(idx,1)/size(m,1)*100;
        varargout{1}=rau(querr,size(m,1));
        meta.ylabel='RAU quadrant errors';
        par.par1=size(idx,1);
        meta.par1='Number of confusions';
        par.par2=size(m,1);
        meta.par2='Number of responses within the lateral range';
      case 'accE'
        varargout{1}=npi2pi(circmean(m(:,4)-m(:,2)));
        meta.ylabel='Elevation bias (deg)';
      case 'absaccE'
        varargout{1}=abs(npi2pi(circmean(m(:,4)-m(:,2))));
        meta.ylabel='Absolute elevation bias (deg)';
      case 'precE'
        precE=circstd(m(:,4)-m(:,2));      
        varargout{1}=precE;
        meta.ylabel='Elevation precision error (deg)';
      case 'lateral'        % obsolete
        nargout=2;
        varargout{1}=accL;
        varargout{2}=precL;
      case 'polar'          % obsolete
        nargout=3;
        varargout{1}=accP;
        varargout{2}=precP;        
        varargout{3}=querr;
      case 'precPmedianlocal'
        idx=find(abs(m(:,7))<=30);
        m=m(idx,:);
        idx=find(abs(npi2pi(m(:,8)-m(:,6)))<90);
        if isempty(idx), error('Quadrant errors of 100%, local bias is NaN'); end
        precPM=circstd(m(idx,8)-m(idx,6));        
        varargout{1}=precPM;
        meta.ylabel='Local central precision error (deg)';
      case 'precPmedian'
        idx=find(abs(m(:,7))<=30);
        m=m(idx,:);
        precPM=circstd(m(:,8)-m(:,6));        
        varargout{1}=precPM;
        meta.ylabel='Central precision error (deg)';
      case 'rmsPmedianlocal'
        idx=find(abs(m(:,7))<=30);
        m=m(idx,:);
        idx=find(abs(npi2pi(m(:,8)-m(:,6)))<90);        
        rmsPlocal=sqrt(sum((npi2pi(m(idx,8)-m(idx,6))).^2)/size(idx,1));
        varargout{1}=rmsPlocal;
        meta.ylabel='Local central RMS error (deg)';
      case 'rmsPmedian'
        idx=find(abs(m(:,7))<=30);
        m=m(idx,:);
        rmsP=sqrt(sum((npi2pi(m(:,8)-m(:,6))).^2)/size(m,1));
        varargout{1}=rmsP;
        meta.ylabel='Central RMS error (deg)';
      case 'SCC'
        trgt = m(:,[1 2]);
        resp = m(:,[3 4]);        
        T = pol2dcos(trgt);
        R = pol2dcos(resp);        
        sxy = T' * R;         
        sxx = T' * T;       
        syy = R' * R;               
        scc = det(sxy) / sqrt(det(sxx) * det(syy));
        varargout{1}=scc;
        meta.ylabel='Spatial correlation coefficient';
      otherwise
        error(['Unknown ERR: ' err]);
    end
    varargout{2}=meta;
    varargout{3}=par; 
  end