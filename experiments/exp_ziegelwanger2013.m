function varargout=exp_ziegelwanger2013(varargin)
%EXP_ZIEGELWANGER2013   Figures from Ziegelwanger and Majdak (2013)
%   Usage: data = exp_ziegelwanger2013(flag)
%
%   `exp_ziegelwanger2013(flags)` reproduces figures of the paper 
%   chapter from Ziegelwanger and Majdak (2013).
%
%   Optional fields of output *data* structure:
%
%   The following flags can be specified:
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'fig2'    Reproduce Fig. 2:
%
%     'fig3'    Reproduce Fig. 3:
%
%     'fig5'    Reproduce Fig. 5:
%
%     'fig6'    Reproduce Fig. 6:
%
%     'fig7'    Reproduce Fig. 7:
%
%     'fig8'    Reproduce Fig. 8:
%
%     'fig9'    Reproduce Fig. 9:
%
%     'fig11'   Reproduce Fig. 11:
%
%     'fig12'   Reproduce Fig. 12:
%
%   See also: ziegelwanger2013, ziegelwanger2013onaxis,
%   ziegelwanger2013offaxis
%
%   Examples:
%   ---------
%
%   To display Fig. 3 use :::
%
%     exp_ziegelwanger2013('fig3');
%
%   References: 

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna, Austria

%% ------ Check input options --------------------------------------------

  definput.flags.type = {'missingflag',...
    'fig2','fig3','fig5','fig6',...
    'fig7','fig8','fig9','fig11','fig12'};
  definput.flags.plot = {'plot','noplot'};

  % Parse input options
  [flags,kv]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% Figure 2
if flags.do_fig2
    
    hpath = which('hrtfinit'); % find local path of hrtf repository
    hpath = hpath(1:end-10);
    
    load([hpath 'ziegelwanger2013/ARI/NH89/hrtf_M_hrtf.mat'])

    subplot(122)
    %---------------------------Threshold---------------------------
    temp=ziegelwanger2013(hM,meta,stimPar,1,0);
    MAX=temp.toa;

    %---------------------------Centroid----------------------------
    temp=ziegelwanger2013(hM,meta,stimPar,2,0);
    CTD=temp.toa;

    %---------------------------Groupdelay--------------------------
    temp=ziegelwanger2013(hM,meta,stimPar,3,0);
    AGD=temp.toa;

    %---------------------------Minimal-Phase-----------------------
    temp=ziegelwanger2013(hM,meta,stimPar,4,0);
    MCM=temp.toa;
    clear temp

    plotziegelwanger2013(MCM(:,1),4,[0 0 0]/255,0,1,1,meta,stimPar,'-',1);
    hold on
    plotziegelwanger2013(MAX(:,1),4,[0 0 80]/255,0,1,1,meta,stimPar,'--',1);
    plotziegelwanger2013(CTD(:,1),4,[50 220 50]/255,0,1,1,meta,stimPar,'-',1);
    plotziegelwanger2013(AGD(:,1),4,[250 80 80]/255,0,1,1,meta,stimPar,'--',1);
    xlim([-10 370])
    ylim([2.65 4.05])
    grid off
    legend(' MCM',' MAX',' CTD',' AGD','Location','NorthWest');
    legend boxoff
    xlabel('Azimuth (°) ')
    ylabel('TOA (ms) ')
%     set(gca,'FontSize',ticklabelsize,'LineWidth',bw,'LineWidth',bw);%,'YTick',[-0.8 -0.4 0 0.4 0.8]
    title('');
    
    subplot(121)
    time=(0:size(hM,1)-1)/stimPar.SamplingRate*1000;
    MAX=round(MAX);
    CTD=round(CTD);
    AGD=round(AGD);
    MCM=round(MCM);
    idx=ARI_FindPosition(meta,85,0);

    fprintf(['MAX: ITD is ' num2str(time(diff(MCM(idx,:),1,2))) ' ms\n'])
    fprintf(['CTD: ITD is ' num2str(time(diff(MAX(idx,:),1,2))) ' ms\n'])
    fprintf(['AGD: ITD is ' num2str(time(diff(CTD(idx,:),1,2))) ' ms\n'])
    fprintf(['MCM: ITD is ' num2str(time(diff(AGD(idx,:),1,2))) ' ms\n'])

    plot(time,hM(:,idx,1)/max(abs(hM(:,idx,1))),'k-');
    hold on
    plot(time,hM(:,idx,2)/max(abs(hM(:,idx,2))),'k--')
    h=stem([time(MCM(idx,1)) time(MCM(idx,1))],[-2 hM(MCM(idx,1),idx,1)/max(abs(hM(:,idx,1)))],'r-','BaseValue',-1);
    stem([time(CTD(idx,1)) time(CTD(idx,1))],[-2 hM(CTD(idx,1),idx,1)/max(abs(hM(:,idx,1)))],'g-','BaseValue',-1)
    stem([time(MAX(idx,1)) time(MAX(idx,1))],[-2 hM(MAX(idx,1),idx,1)/max(abs(hM(:,idx,1)))],'b-','BaseValue',-1)
    stem([time(AGD(idx,1)) time(AGD(idx,1))],[-2 hM(AGD(idx,1),idx,1)/max(abs(hM(:,idx,1)))],'m-','BaseValue',-1)
    stem([time(MCM(idx,2)) time(MCM(idx,2))],[-2 hM(MCM(idx,2),idx,2)/max(abs(hM(:,idx,2)))],'r-','BaseValue',-1)
    stem([time(CTD(idx,2)) time(CTD(idx,2))],[-2 hM(CTD(idx,2),idx,2)/max(abs(hM(:,idx,2)))],'g-','BaseValue',-1)
    stem([time(AGD(idx,2)) time(AGD(idx,2))],[-2 hM(AGD(idx,2),idx,2)/max(abs(hM(:,idx,2)))],'m-','BaseValue',-1)
    stem([time(MAX(idx,2)) time(MAX(idx,2))],[-2 hM(MAX(idx,2),idx,2)/max(abs(hM(:,idx,2)))],'b-','BaseValue',-1)
    plot(time(MAX(idx,1)),hM(MAX(idx,1),idx,1)/max(abs(hM(:,idx,1))),'b^','MarkerFaceColor','b')
    plot(time(MCM(idx,1)),hM(MCM(idx,1),idx,1)/max(abs(hM(:,idx,1))),'ro','MarkerFaceColor','r')
    plot(time(CTD(idx,1)),hM(CTD(idx,1),idx,1)/max(abs(hM(:,idx,1))),'gd','MarkerFaceColor','g')
    plot(time(AGD(idx,1)),hM(AGD(idx,1),idx,1)/max(abs(hM(:,idx,1))),'ms','MarkerFaceColor','m')
    plot(time(MAX(idx,2)),hM(MAX(idx,2),idx,2)/max(abs(hM(:,idx,2))),'b^','MarkerFaceColor','b')
    plot(time(MCM(idx,2)),hM(MCM(idx,2),idx,2)/max(abs(hM(:,idx,2))),'ro','MarkerFaceColor','r')
    plot(time(CTD(idx,2)),hM(CTD(idx,2),idx,2)/max(abs(hM(:,idx,2))),'gd','MarkerFaceColor','g')
    plot(time(AGD(idx,2)),hM(AGD(idx,2),idx,2)/max(abs(hM(:,idx,2))),'ms','MarkerFaceColor','m')
    xlim([2.4 4.1])
    ylim([-1.1 1.1])
    xlabel('Time (ms) ')
    ylabel('Amplitude ')
%     set(gca,'FontSize',ticklabelsize,'LineWidth',bw,'LineWidth',bw,'XTick',[2.5 3 3.5 4]);
    set(get(h,'Baseline'),'Visible','off')
end

%% Figure 3, Figure 6
if flags.do_fig3 || flags.do_fig6
    
    hpath = which('hrtfinit'); % find local path of hrtf repository
    hpath = hpath(1:end-10);
    
    load([hpath 'ziegelwanger2013/ARI/NH89/hrtf_M_hrtf.mat'])

    p0_onaxis=[[0.0875; pi/2; 0; 0.0001] [0.0875; -pi/2; 0; 0.0001]];
    p0_onaxis=transpose(p0_onaxis);
    p_onaxis=zeros(size(p0_onaxis));
    p0_offaxis=zeros(2,7);
    p_offaxis=p0_offaxis;
    global ch

    toa=zeros(size(hM,2),size(hM,3));
    toaEst=zeros(size(hM,2),size(hM,3));
    indicator1=zeros(size(hM,2),size(hM,3));
    indicator3=indicator1;
    meta.pos(:,8)=cumsum(ones(size(meta.pos,1),1));
    hM_min=ARI_MinimalPhase([hM; zeros(4096-size(hM,1),size(hM,2),size(hM,3))]);
    hM_min=hM_min(1:size(hM,1),:,:);
    for ii=1:size(hM,2)
        for jj=1:size(hM,3)
            if isnan(hM_min(1,ii,jj))
                hM_min(:,ii,jj)=ARI_MinimalPhase(hM(:,ii,jj));
            end
        end
    end
    corrcoeff=zeros(size(hM,2),size(hM,3));
    for ii=1:size(hM,2)
        for jj=1:size(hM,3)
            [c,lag]=xcorr(transpose(squeeze(hM(:,ii,jj))),transpose(squeeze(hM_min(:,ii,jj))),size(hM,1)-1,'none');
            [corrcoeff(ii,jj),idx]=max(abs(c));
            corrcoeff(ii,jj)=corrcoeff(ii,jj)/sum(hM(:,ii,jj).^2);
            toaEst(ii,jj)=lag(idx);
        end
    end
    
    for ch=1:size(hM,3)
        % Outlier detection: smooth TOA in horizontal planes
        epsilon=5;
        slope=zeros(size(hM,2),1);
        for ele=min(meta.pos(:,2)):epsilon:max(meta.pos(:,2)) %calculate slope for each elevation along azimut
            idx=find(meta.pos(:,2)>ele-epsilon/2 & meta.pos(:,2)<=ele+epsilon/2);
            if numel(idx)>1
                idx(length(idx)+1)=idx(1);
                slope(idx(1:end-1),1)=diff(toaEst(idx,ch))./abs(diff(meta.pos(idx,1)));
            end
        end
        sloperms=sqrt(sum(slope.^2)/length(slope));
        if sloperms<30/(length(find(meta.pos(:,2)==0))/2)
            sloperms=30/(length(find(meta.pos(:,2)==0))/2);
        end
        for ele=min(meta.pos(:,2)):epsilon:max(meta.pos(:,2))
            idx=find(meta.pos(:,2)>ele-epsilon/2 & meta.pos(:,2)<=ele+epsilon/2);
            for ii=1:length(idx)-1
                if abs(slope(idx(ii)))>sloperms
                    for jj=0:1
                        if ii+jj==0 || ii+jj==length(idx)
                            indicator1(idx(end),ch)=1;
                        else
                            indicator1(idx(mod(ii+jj,length(idx))),ch)=1;
                        end
                    end
                end
            end
            clear idx
        end
        indicator2=indicator1(:,ch);

        % Outlier detection: constant TOA in sagittal planes
        epsilon=2;
        for ii=1:20
            lin=zeros(size(hM,2),1);
            for lat=-90:epsilon:90
                idx=find(meta.pos(:,6)>lat-epsilon/2 & meta.pos(:,6)<=lat+epsilon/2 & indicator2==0);
                idx2=find(meta.pos(:,6)>lat-epsilon/2 & meta.pos(:,6)<=lat+epsilon/2 & indicator1(:,ch)==0);
                if length(idx2)>2
                    lin(idx,1)=toaEst(idx,ch)-mean(toaEst(idx2,ch));
                end
            end
            linrms=sqrt(sum(lin.^2)/length(lin));
            if linrms<2
                linrms=2;
            end
            indicator1(:,ch)=zeros(size(indicator1,1),1);
            indicator3(:,ch)=zeros(size(indicator3,1),1);
            indicator3(abs(lin)>linrms,ch)=ones(length(find(abs(lin)>linrms)),1);
            indicator1(abs(lin)>linrms | indicator2==1,ch)=ones(length(find(abs(lin)>linrms | indicator2==1)),1);
        end

        if flags.do_fig3 && ch==1 %Figure 3
            subplot(121)
            plot([-90; 270],[linrms/stimPar.SamplingRate*1000000; linrms/stimPar.SamplingRate*1000000],'r--');
            hold on
            plot(real(meta.pos(:,7)),abs(lin)/stimPar.SamplingRate*1000000,'b.');
            xlim([-98 278])
            ylim([-5 max(abs(lin)/stimPar.SamplingRate*1000000)+5])
            xlabel('Polar angle (°)')
            ylabel('Sagittal TOA deviation (µs)')
            title('')
            set(gca,'XTick',[-90 0 90 180 270])
        end
        clear lin; clear linrms;

        if flags.do_fig3 && ch==1 %Figure 3
            subplot(122)
            plotziegelwanger2013(toaEst,4,'k',0,1,1,meta,stimPar,{'-'},1);
            hold on
            h=plotziegelwanger2013(indicator3.*toaEst,3,'w',0,1,1,meta,stimPar,{'^'},4);
            set(h,'MarkerFaceColor','w','MarkerEdgeColor','w');
            h=plotziegelwanger2013(indicator2*[1 1].*toaEst,3,'w',0,1,1,meta,stimPar,{'v'},4);
            set(h,'MarkerFaceColor','w','MarkerEdgeColor','w');
            h=plotziegelwanger2013(indicator3.*toaEst,3,'b',0,1,1,meta,stimPar,{'^'},4);
            set(h,'LineWidth',2);
            h=plotziegelwanger2013(indicator2*[1 1].*toaEst,3,'b',0,1,1,meta,stimPar,{'v'},4);
            set(h,'LineWidth',2);
            h=plotziegelwanger2013((-indicator1+1).*toaEst,3,'r',0,1,1,meta,stimPar,{'o'},4);
            set(h,'MarkerFaceColor','r','MarkerEdgeColor','r');
            ylabel('TOA (ms) ')
            xlabel('Azimuth (°) ')
            grid off
            xlim([-10 370])
            ylim([2.65 3.65])
            title('')
            set(gca,'YTick',[2.7 3 3.3 3.6])
        end
    end

    for ch=1:size(hM,3)
        p0_onaxis(ch,4)=min(toaEst(indicator1(:,ch)==0,ch))/stimPar.SamplingRate;
        p0off_onaxis=[0.06 pi/4 pi/4 0.001];

        % Fit on-axis model to outlier adjusted set of estimated TOAs
        idx=find(indicator1(:,ch)==0);
        x=meta.pos(idx,1:2)*pi/180;
        y=toaEst(idx,ch)/stimPar.SamplingRate;
        p_onaxis(ch,:)=lsqcurvefit(@ziegelwanger2013onaxis,p0_onaxis(ch,:),x,y,p0_onaxis(ch,:)-p0off_onaxis,p0_onaxis(ch,:)+p0off_onaxis,optimset('Display','off','TolFun',1e-6));%1e-4 for spheres
        toa(:,ch)=ziegelwanger2013onaxis(p_onaxis(ch,:),meta.pos(:,1:2)*pi/180)*stimPar.SamplingRate;
    end

    TolFun=[1e-5 1e-6; 1e-6 1e-8];
    % Fit off-axis model to outlier adjusted set of estimated TOAs
    for ii=1:size(TolFun,1)
        for ch=1:size(hM,3)
            idx=find(indicator1(:,ch)==0);
            x=meta.pos(idx,1:2)*pi/180;
            y=toaEst(idx,ch)/stimPar.SamplingRate;
            p0_offaxis(ch,:)=[p0_onaxis(ch,1) 0 0 0 p0_onaxis(ch,4) p0_onaxis(ch,2) p0_onaxis(ch,3)];
            p0off_offaxis=[0.05 0.05 0.05 0.05 0.001 pi pi];
            p_offaxis(ch,:)=lsqcurvefit(@ziegelwanger2013offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0off_offaxis,p0_offaxis(ch,:)+p0off_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));%1e-6 for spheres
            toa(:,ch)=ziegelwanger2013offaxis(p_offaxis(ch,:),meta.pos(:,1:2)*pi/180)*stimPar.SamplingRate;
        end
        if abs(diff(p_offaxis(:,1)))>0.003 || abs(diff(p_offaxis(:,3)))>0.003
            p_offaxis(:,[1 3])=p_offaxis([2 1],[1 3]);
            for ch=1:size(hM,3)
                idx=find(indicator1(:,ch)==0);
                x=meta.pos(idx,1:2)*pi/180;
                y=toaEst(idx,ch)/stimPar.SamplingRate;
                p0_offaxis(ch,:)=[p_offaxis(ch,1) mean(p_offaxis(:,2)) p_offaxis(ch,3) mean(p_offaxis(:,4)) mean(p_offaxis(:,5)) p_offaxis(ch,6) p_offaxis(ch,7)];
                p0off_offaxis=[0.05 0.05 0.05 0.05 0.001 pi/2 pi/2];
                p_offaxis(ch,:)=lsqcurvefit(@ziegelwanger2013offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0off_offaxis,p0_offaxis(ch,:)+p0off_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));%1e-8 for spheres 
                toa(:,ch)=ziegelwanger2013offaxis(p_offaxis(ch,:),meta.pos(:,1:2)*pi/180)*stimPar.SamplingRate;
            end
        end
        if abs(diff(p_offaxis(:,1)))<0.003 && abs(diff(p_offaxis(:,2)))<0.003 && abs(diff(p_offaxis(:,3)))<0.003 && abs(diff(p_offaxis(:,4)))<0.003
            break
        end
    end

    if flags.do_fig6 %Figure 6
        h=plotziegelwanger2013(indicator1.*toaEst,3,'b',0,1,1,meta,stimPar,{'^'},4);
        set(h,'LineWidth',2);
        hold on
        h=plotziegelwanger2013(indicator1.*toaEst,3,'b',0,1,1,meta,stimPar,{'v'},4);
        set(h,'LineWidth',2);
        h=plotziegelwanger2013((-indicator1+1).*toaEst,3,'r',0,1,1,meta,stimPar,{'o'},4);
        set(h,'MarkerFaceColor','r','MarkerEdgeColor','r');
        plotziegelwanger2013(toa,4,'k',0,1,1,meta,stimPar,{'-'},1);
        ylabel('TOA (ms) ')
        xlabel('Azimuth (°) ')
        grid off
        xlim([-10 370])
        ylim([2.65 3.65])
        title('')
        set(gca,'YTick',[2.7 3 3.3 3.6])
    end
end

%% Figure 5
if flags.do_fig5
    
    hpath = which('hrtfinit'); % find local path of hrtf repository
    hpath = hpath(1:end-10);
    
    s=load([hpath 'ziegelwanger2013' filesep 'Sphere' filesep 'Rotation' filesep 'subjects.mat']);
    p=zeros(4,2,length(s.phi));
    
    for ii=1:length(s.subjects)
        load([hpath 'ziegelwanger2013' filesep 'Sphere' filesep 'Rotation' filesep s.subjects{ii} filesep 'hrtf_M_hrtf.mat']);
        meta=ziegelwanger2013(hM,meta,stimPar,4,1);
        p(:,:,ii)=meta.p_onaxis;
    end

    sym='sdo'; %plot symbols
    clr=[0,0,255; 255,0,0; 255,255,67]/255; %plot colors
    meclr=[0,0,255; 255,0,0; 255,255,67]/255; %marker edge colors
    
    % radii
    subplot(311)
    var=[squeeze(p(1,1,:))*100 squeeze(p(1,2,:))*100 s.radius(:)/10];
    err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
    err=reshape(err,numel(err),1);
    fprintf(['Radius: average err is ' num2str(mean(err)) ' cm\n']);
    fprintf(['        standard deviation is ' num2str(std(err)) ' cm\n']);
    for ch=1:size(p,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    plot(var(:,3),'k--')
    clear var;
    ylabel('r (cm) ')

    %phi
    subplot(312)
    var=[squeeze(p(2,1,:))/pi*180 squeeze(p(2,2,:))/pi*180 -s.phi+ones(length(s.phi),1)*90 -s.phi-ones(length(s.phi),1)*90];
    err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,4)]);
    err=reshape(err,numel(err),1);
    fprintf(['Phi: average err is ' num2str(mean(err)) '°\n']);
    fprintf(['     standard deviation is ' num2str(std(err)) '°\n']);
    for ch=1:size(p,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    for ch=1:size(p,2)
        plot(1:9,var(1:9,2+ch),'k--')
        plot(10:14,var(10:14,2+ch),'k--')
        plot(15:23,var(15:23,2+ch),'k--')
        plot(24:28,var(24:28,2+ch),'k--')
        plot(29:37,var(29:37,2+ch),'k--')
        plot(38:42,var(38:42,2+ch),'k--')
    end
    clear var;
    set(gca,'YTick',[-90 90])
    ylabel('\phi_e (°) ')

    %theta
    subplot(313)
    var=[squeeze(p(3,1,:))/pi*180 squeeze(p(3,2,:))/pi*180 s.theta -s.theta];
    err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,4)]);
    err=reshape(err,numel(err),1);
    fprintf(['Theta: average err is ' num2str(mean(err)) '°\n']);
    fprintf(['       standard deviation is ' num2str(std(err)) '°\n']);
    for ch=1:size(p,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    for ch=1:size(p,2)
        plot(1:9,var(1:9,2+ch),'k--')
        plot(10:14,var(10:14,2+ch),'k--')
        plot(15:23,var(15:23,2+ch),'k--')
        plot(24:28,var(24:28,2+ch),'k--')
        plot(29:37,var(29:37,2+ch),'k--')
        plot(38:42,var(38:42,2+ch),'k--')
    end
    clear var;
    ylabel('\theta_e (°) ')
    xlabel('Condition')
end

%% Figure 7, Figure 12
if flags.do_fig7 || flags.do_fig12
    hrtf={'ARI','CIPIC','LISTEN'};

    sym='ods'; %plot symbols
    
    hpath = which('hrtfinit'); % find local path of hrtf repository
    hpath = hpath(1:end-10);

    %-------------------------------Load Data----------------------------------
    for kk=1:length(hrtf)
        load(directoryCheck([hpath 'ziegelwanger2013' filesep hrtf{kk} filesep 'results.mat']))
        if kk==3
            results=results([1:27 29:end]);
        end
        temp1=zeros(size(results(1).meta.p,1),size(results(1).meta.p,2),length(results));
        temp2=zeros(size(results(1).meta.pExt,1),size(results(1).meta.pExt,2),length(results));
        temp3=zeros(size(results(1).meta.pExtFull,1),size(results(1).meta.pExtFull,2),length(results));
        temp4=zeros(length(results),4);
        temp6=zeros(length(results),4);
        temp7=zeros(length(results),4);
        for ii=1:length(results)
            temp1(:,:,ii)=results(ii).meta.p;
            temp2(:,1:size(results(ii).meta.pExt,2),ii)=results(ii).meta.pExt;
            temp3(:,1:size(results(ii).meta.pExtFull,2),ii)=results(ii).meta.pExtFull;
            temp4(ii,:)=[results(ii).meta.performance(1).outliers results(ii).meta.performance(2).outliers results(ii).meta.performance(3).outliers results(ii).meta.performance(4).outliers];
            if isfield(results(ii).meta.performance(1)','outliersl')
                temp6(ii,:)=[results(ii).meta.performance(1).outliersl results(ii).meta.performance(2).outliersl results(ii).meta.performance(3).outliersl results(ii).meta.performance(4).outliersl];
                temp7(ii,:)=[results(ii).meta.performance(1).outliersr results(ii).meta.performance(2).outliersr results(ii).meta.performance(3).outliersr results(ii).meta.performance(4).outliersr];
            end
        end
        p_onaxis{kk}=temp1;
        p_offaxis{kk}=temp3;
        outliers{kk}=temp4;
        outliersear{kk}=[temp6; temp7];
    end
    
    if flags.do_fig7 %Figure 7
        fprintf('On-Axis Model:\n')
        fprintf(['Average Radius: ' num2str(mean([mean(squeeze(p_onaxis{1}(1,:,:)*100)) mean(squeeze(p_onaxis{2}(1,:,:)*100)) mean(squeeze(p_onaxis{3}(1,:,:)*100))])) ' cm\n'])
        fprintf(['Average Radius Difference: ' num2str(mean([mean(abs(diff(squeeze(p_onaxis{1}(1,:,:)*100)))) mean(abs(diff(squeeze(p_onaxis{2}(1,:,:)*100)))) mean(abs(diff(squeeze(p_onaxis{3}(1,:,:)*100))))])) 'cm\n'])
        fprintf(['Maximum Radius Difference: ' num2str(max([max(abs(diff(squeeze(p_onaxis{1}(1,:,:)*100)))) max(abs(diff(squeeze(p_onaxis{2}(1,:,:)*100)))) max(abs(diff(squeeze(p_onaxis{3}(1,:,:)*100))))])) ' cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(mean(p_onaxis{kk}(1,:,:)*100,2)) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_onaxis{kk}(1,1,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_onaxis{kk}(1,2,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            temp=size(var,1)+1;
        end
        [~,idx]=sort(var(:,1));
        var=var(idx,:);
        varl=varl(idx,:);
        varr=varr(idx,:);
        var(:,3)=transpose(1:size(var,1));
        varl(:,3)=transpose(1:size(var,1));
        varr(:,3)=transpose(1:size(var,1));
        for ii=1:size(varl,1)
            stm=stem(ii,varl(ii,1),'--k','BaseValue',varr(ii,1));
            hold on
        end
        baseline_handle = get(stm,'BaseLine');
        set(baseline_handle,'LineStyle','none')
        for kk=1:length(hrtf)
            h{kk}=plot(varl(varl(:,2)==kk,3),varl(varl(:,2)==kk,1),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
            plot(varr(varr(:,2)==kk,3),varr(varr(:,2)==kk,1),sym(kk),'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        clear var varl varr stm;
        xlabel('Listeners ')
        ylabel('r (cm) ')
        xlim([-1.5 temp+1.5])
        ylim([4.1 12+0.7])
        legend([h{1},h{2},h{3}],'ARI','CIPIC','LISTEN','Location','NorthWest');
        legend boxoff
    end
    
    if flags.do_fig12 %Figure12
        %radii
        subplot(411)
        fprintf('Off-Axis Model:\n')
        fprintf(['Average Radius: ' num2str(mean([mean(squeeze(p_offaxis{1}(1,:,:)*100)) mean(squeeze(p_offaxis{2}(1,:,:)*100)) mean(squeeze(p_offaxis{3}(1,:,:)*100))])) ' cm\n'])
        fprintf(['Average Radius Difference: ' num2str(mean([mean(abs(diff(squeeze(p_offaxis{1}(1,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{2}(1,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{3}(1,:,:)*100))))])) 'cm\n'])
        fprintf(['Maximum Radius Difference: ' num2str(max([max(abs(diff(squeeze(p_offaxis{1}(1,:,:)*100)))) max(abs(diff(squeeze(p_offaxis{2}(1,:,:)*100)))) max(abs(diff(squeeze(p_offaxis{3}(1,:,:)*100))))])) ' cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(mean(p_offaxis{kk}(1,:,:)*100,2)) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(1,1,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(1,2,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            temp=size(var,1)+1;
        end
        [~,idx]=sort(var(:,1));
        var=var(idx,:);
        varl=varl(idx,:);
        varr=varr(idx,:);
        var(:,3)=transpose(1:size(var,1));
        varl(:,3)=transpose(1:size(var,1));
        varr(:,3)=transpose(1:size(var,1));
        for ii=1:size(varl,1)
            stm=stem(ii,varl(ii,1),'--k','BaseValue',varr(ii,1));
            hold on
        end
        baseline_handle = get(stm,'BaseLine');
        set(baseline_handle,'LineStyle','none')
        for kk=1:length(hrtf)
            h{kk}=plot(varl(varl(:,2)==kk,3),varl(varl(:,2)==kk,1),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
            plot(varr(varr(:,2)==kk,3),varr(varr(:,2)==kk,1),sym(kk),'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        clear var varl varr;
        ylabel('r (cm) ')
        xlim([-1.5 temp+1.5])
        ylim([4.1 12+0.7])
        legend([h{1},h{2},h{3}],'ARI','CIPIC','LISTEN','Location','NorthWest');
        legend boxoff

        %lateral displacements
        % xM
        subplot(412)
        fprintf(['Average xM: ' num2str(mean([mean(abs(squeeze(p_offaxis{1}(2,:,:)*100))) mean(abs(squeeze(p_offaxis{2}(2,:,:)*100))) mean(abs(squeeze(p_offaxis{3}(2,:,:)*100)))])) ' cm\n'])
        fprintf(['Average xM Difference: ' num2str(mean([mean(abs(diff(squeeze(p_offaxis{1}(2,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{2}(2,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{3}(2,:,:)*100))))])) 'cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(mean(p_offaxis{kk}(2,:,:)*100,2)) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(2,1,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(2,2,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            temp=size(var,1)+1;
        end
        var=var(idx,:);
        varl=varl(idx,:);
        varr=varr(idx,:);
        var(:,3)=transpose(1:size(var,1));
        varl(:,3)=transpose(1:size(var,1));
        varr(:,3)=transpose(1:size(var,1));
        for ii=1:size(varl,1)
            stm=stem(ii,varl(ii,1),'--k','BaseValue',varr(ii,1));
            hold on
        end
        baseline_handle = get(stm,'BaseLine');
        set(baseline_handle,'LineStyle','none')
        for kk=1:length(hrtf)
            h{kk}=plot(varl(varl(:,2)==kk,3),varl(varl(:,2)==kk,1),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
            plot(varr(varr(:,2)==kk,3),varr(varr(:,2)==kk,1),sym(kk),'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        clear var varl varr;
        ylabel('x_M (cm) ')
        xlim([-1.5 temp+1.5])
        ylim([-4.5 4.5])

        % yM
        subplot(413)
        fprintf(['Average yM: ' num2str(mean([mean(abs(squeeze(p_offaxis{1}(3,:,:)*100))) mean(abs(squeeze(p_offaxis{2}(3,:,:)*100))) mean(abs(squeeze(p_offaxis{3}(3,:,:)*100)))])) ' cm\n'])
        fprintf(['Average yM Difference: ' num2str(mean([mean(abs(diff(squeeze(p_offaxis{1}(3,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{2}(3,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{3}(3,:,:)*100))))])) 'cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(mean(p_offaxis{kk}(3,:,:)*100,2)) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(3,1,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(3,2,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            temp=size(var,1)+1;
        end
        var=var(idx,:);
        varl=varl(idx,:);
        varr=varr(idx,:);
        var(:,3)=transpose(1:size(var,1));
        varl(:,3)=transpose(1:size(var,1));
        varr(:,3)=transpose(1:size(var,1));
        for ii=1:size(varl,1)
            stm=stem(ii,varl(ii,1),'--k','BaseValue',varr(ii,1));
            hold on
        end
        baseline_handle = get(stm,'BaseLine');
        set(baseline_handle,'LineStyle','none')
        for kk=1:length(hrtf)
            h{kk}=plot(varl(varl(:,2)==kk,3),varl(varl(:,2)==kk,1),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
            plot(varr(varr(:,2)==kk,3),varr(varr(:,2)==kk,1),sym(kk),'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        clear var;
        ylabel('y_M (cm) ')
        xlim([-1.5 temp+1.5])
        ylim([-3.5 4.5])

        % zM
        subplot(414)
        fprintf(['Average zM: ' num2str(mean([mean(abs(squeeze(p_offaxis{1}(4,:,:)*100))) mean(abs(squeeze(p_offaxis{2}(4,:,:)*100))) mean(abs(squeeze(p_offaxis{3}(4,:,:)*100)))])) ' cm\n'])
        fprintf(['Average zM Difference: ' num2str(mean([mean(abs(diff(squeeze(p_offaxis{1}(4,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{2}(4,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{3}(4,:,:)*100))))])) 'cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(mean(p_offaxis{kk}(4,:,:)*100,2)) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(4,1,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(4,2,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
            temp=size(var,1)+1;
        end
        var=var(idx,:);
        varl=varl(idx,:);
        varr=varr(idx,:);
        var(:,3)=transpose(1:size(var,1));
        varl(:,3)=transpose(1:size(var,1));
        varr(:,3)=transpose(1:size(var,1));
        for ii=1:size(varl,1)
            stm=stem(ii,varl(ii,1),'--k','BaseValue',varr(ii,1));
            hold on
        end
        baseline_handle = get(stm,'BaseLine');
        set(baseline_handle,'LineStyle','none')
        for kk=1:length(hrtf)
            h{kk}=plot(varl(varl(:,2)==kk,3),varl(varl(:,2)==kk,1),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
            plot(varr(varr(:,2)==kk,3),varr(varr(:,2)==kk,1),sym(kk),'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        clear var;
        xlabel('Listeners')
        ylabel('z_M (cm) ')
        xlim([-1.5 temp+1.5])
        ylim([-4.5 6.5])
    end
end

%% Figure 8
if flags.do_fig8
    
    hpath = which('hrtfinit'); % find local path of hrtf repository
    hpath = hpath(1:end-10);
    
    load([hpath 'ziegelwanger2013/ARI/NH89/hrtf_M_hrtf.mat'])
    
    temp=ziegelwanger2013(hM,meta,stimPar,4,1);
    meta.toa5=temp.toa;
    fprintf(['Radii for left and right ear: ' num2str(temp.p_onaxis(1,:)*100) ' cm\n'])
    temp=TOA_Calc(hM,meta,stimPar,4,0,0);
    meta.toa6=temp.toa;
    fprintf(['Maximum TOA difference left: ' num2str((max(meta.toa5(:,1))-min(meta.toa5(:,1)))/stimPar.SamplingRate*1000) ' ms\n'])
    fprintf(['Maximum TOA difference right: ' num2str((max(meta.toa5(:,2))-min(meta.toa5(:,2)))/stimPar.SamplingRate*1000) ' ms\n'])

    plotziegelwanger2013(meta.toa6-min(min(meta.toa5)),4,'k',0,1,1,meta,stimPar,'--',1);
    plotziegelwanger2013(meta.toa5-min(min(meta.toa5)),4,'k',0,1,1,meta,stimPar,'--',2);
    plotziegelwanger2013(meta.toa6-min(min(meta.toa5)),4,'k',0,2,1,meta,stimPar,'-',1);
    plotziegelwanger2013(meta.toa5-min(min(meta.toa5)),4,'k',0,2,1,meta,stimPar,'-',2);

    xlim([-5 365])
    ylim([-0.05 0.95])
    grid off
    xlabel('Azimuth (°) ')
    ylabel('Relative TOA (ms) ');
    title('')
end

%% Figure 9, Figure 11
if flags.do_fig9 || flags.do_fig11
    sym='sdo'; %plot symbols
    clr=[0,0,255; 255,0,0; 255,255,67]/255; %plot colors
    meclr=[0,0,255; 255,0,0; 255,255,67]/255; %marker edge colors
    
    hpath = which('hrtfinit'); % find local path of hrtf repository
    hpath = hpath(1:end-10);
    
    s=load([hpath 'ziegelwanger2013' filesep 'Sphere' filesep 'Displacement' filesep 'subjects.mat']);
    
    p_onaxis=zeros(4,2,length(s.subjects));
    p_offaxis=zeros(7,2,length(s.subjects));
    for ii=1:length(s.subjects)
        load([hpath 'ziegelwanger2013' filesep 'Sphere' filesep 'Displacement' filesep s.subjects{ii} filesep 'hrtf_M_hrtf.mat']);
        meta=ziegelwanger2013(hM,meta,stimPar,4,1);
        p_onaxis(:,:,ii)=meta.p_onaxis;
        p_offaxis(:,:,ii)=meta.p_offaxis;
    end
    
    if flags.do_fig9 %Figure 9
        p1=p_onaxis(:,:,[1:3 length(s.xM)/3+1:length(s.xM)/3+3 length(s.xM)/3*2+1:length(s.xM)/3*2+3]);
        r1=s.radius([1:3 length(s.xM)/3+1:length(s.xM)/3+3 length(s.xM)/3*2+1:length(s.xM)/3*2+3]);
        yM1=s.yM([1:3 length(s.xM)/3+1:length(s.xM)/3+3 length(s.xM)/3*2+1:length(s.xM)/3*2+3]);
        
        %radii
        subplot(211)
%         ymax2=round(max(max(squeeze(p1(1,:,:)*100))))+1;
        var=[squeeze(p1(1,1,:))*100 squeeze(p1(1,2,:))*100 r1/10];
        for ch=1:size(p1,2)
            plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
            hold on
        end
        plot(r1/10,'k--')
        clear var;
        ylabel('r (cm) ')

        %yM
        subplot(212)
        plot(-yM1*100,'k--')
        clear var;
        xlabel('Condition')
        ylabel('y_M (cm) ')
    end
    
    if flags.do_fig11 %Figure 11
        center=[s.xM(1:length(s.xM)/3) s.yM(1:length(s.yM)/3) s.zM(1:length(s.zM)/3)];
        [~,idx3]=sort(squeeze(center(3,:,1)));
        [~,idx2]=sort(squeeze(center(1,:,1)));
        [~,idx1]=sort(squeeze(center(2,:,1)));
        idx=idx3(idx2(idx1));
        idx=[idx; idx+length(xM)/3; idx+length(xM)/3*2];
        s.radius=s.radius(idx);
        p_offaxis=p_offaxis(:,:,idx);
        s.xM=s.xM(idx);
        s.yM=s.yM(idx);
        s.zM=s.zM(idx);

        %radii
        subplot(411)
        temp=1;
        ymax2=round(max(max(squeeze(p_offaxis(1,:,:)*100))))+1;
        var=[squeeze(p_offaxis(1,1,:))*100 squeeze(p_offaxis(1,2,:))*100 s.radius/10];
        err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
        err=reshape(err,numel(err),1);
        temp1=err;
        fprintf(['Radius: average err is ' num2str(mean(err)) ' cm\n']);
        fprintf(['        standard deviation is ' num2str(std(err)) ' cm\n']);
        fprintf(['        maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
        fprintf(['        average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        for ch=1:size(pExtp_offaxisFull,2)
            plot(1:temp-1,var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
            hold on
        end
        plot(var(:,3),'k--')
        clear var;
        ylabel('r (cm) ')
        xlim([-1 temp+1])
        ylim([4.5 ymax2])

        %xM
        subplot(412)
        var=[squeeze(p_offaxis(2,1,:))*100 squeeze(p_offaxis(2,2,:))*100 -s.xM*100];
        err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
        err=reshape(err,numel(err),1);
        fprintf(['xM: average err is ' num2str(mean(err)) ' cm\n']);
        fprintf(['    standard deviation is ' num2str(std(err)) ' cm\n']);
        fprintf(['    maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
        fprintf(['    average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        for ch=1:size(p_offaxis,2)
            plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
            hold on
        end  
        plot(var(:,3),'k--')
        clear var;
        ylabel('x_M (cm) ')
        xlabel([-1 temp+1])
        ylabel([-2.5 0.5])

        %yM
        subplot(413)
        var=[squeeze(p_offaxis(3,1,:))*100 squeeze(p_offaxis(3,2,:))*100 -s.yM*100];
        err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
        err=reshape(err,numel(err),1);
        temp1=[temp1; err];
        fprintf(['yM: average err is ' num2str(mean(err)) ' cm\n']);
        fprintf(['    standard deviation is ' num2str(std(err)) ' cm\n']);
        fprintf(['    maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
        fprintf(['    average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        for ch=1:size(p_offaxis,2)
            plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
            hold on
        end
        plot(var(:,3),'k--')
        clear var;
        ylabel('y_M (cm) ')
        xlim([-1 temp+1])
        ylim([-0.5 2.5])

        %zM
        subplot(414)
        var=[squeeze(p_offaxis(4,1,:))*100 squeeze(p_offaxis(4,2,:))*100 -s.zM*100];
        err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
        err=reshape(err,numel(err),1);
        temp1=[temp1; err];
        fprintf(['zM: average err is ' num2str(mean(err)) ' cm\n']);
        fprintf(['    standard deviation is ' num2str(std(err)) ' cm\n']);
        fprintf(['    maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
        fprintf(['    average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        for ch=1:size(p_offaxis,2)
            plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
            hold on
        end
        plot(var([1 3],3),'k--')
        plot(var([2 4],3),'k--')
        clear var;
        set(gca,'YTick',[-1 0])
        xlabel('Condition')
        ylabel('z_M (cm) ')
        xlim([-1 temp+1])
        ylim([-1.5 0.5])

        
        fprintf(['offset: average err is ' num2str(mean(temp1)) ' cm\n']);
        fprintf(['        standard deviation is ' num2str(std(temp1)) ' cm\n']);
        fprintf(['        maximum is ' num2str(max(abs(temp1))) ' cm\n']);
    end
end

end

function idx=ARI_FindPosition(meta,azimuth,elevation)
    psi=sin(elevation/180*pi).*sin(meta.pos(:,2)/180*pi) + ...
        cos(elevation/180*pi).*cos(meta.pos(:,2)/180*pi).*...
        cos(azimuth/180*pi-meta.pos(:,1)/180*pi);
    [~,idx]=min(acos(psi));
end

function out=ARI_MinimalPhase(in)
    n=size(in,1);
    itnr=size(in,2);
    rec=size(in,3);
    out=zeros(size(in));

    for jj=1:rec
        for ii=1:itnr
            h=squeeze(in(:,ii,jj));
            % decompose signal
            amp1=abs(fft(h));

            % transform
            amp2=amp1;
            an2u=-imag(hilbert(log(amp1))); % minimal phase

            % reconstruct signal from amp2 and an2u
            % build a symmetrical phase 
            an2u=an2u(1:floor(n/2)+1);
            an2u=[an2u; -flipud(an2u(2:end+mod(n,2)-1))];
            an2=an2u-round(an2u/2/pi)*2*pi;  % wrap around +/-pi: wrap(x)=x-round(x/2/pi)*2*pi
            % amplitude
            amp2=amp2(1:floor(n/2)+1);
            amp2=[amp2; flipud(amp2(2:end+mod(n,2)-1))];
            % back to time domain
            h2=real(ifft(amp2.*exp(1i*an2)));
            out(:,ii,jj)=h2;
        end
    end
end