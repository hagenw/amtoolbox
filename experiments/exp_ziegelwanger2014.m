function varargout=exp_ziegelwanger2014(varargin)
%EXP_ZIEGELWANGER2013   Figures from Ziegelwanger and Majdak (2013)
%   Usage: data = exp_ziegelwanger2014(flag)
%
%   `exp_ziegelwanger2014(flags)` reproduces figures of the paper from
%   Ziegelwanger and Majdak (2013).
%
%   Optional fields of output *data* structure:
%
%   The following flags can be specified:
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, print results in the console.
% 
%     'reload'	Reload previously calculated results. This is the default.
%
%     'recalc'	Recalculate results.
%
%     'fig2'    Reproduce Fig. 2:
%               
%               Left panel: 
%               Normalized HRIRs of NH89 (ARI database). Sound source was
%               placed 45Â° left in the horizontal plane.
%               Solid line: for the left ear
%               Dashed line: for the right ear
%               Vertical lines: Estimated TOAs
%               Arrows: Resulting ITDs from the corresponding TOAs
%               Circles, Red: Minimum-Phase Cross-Correlation Method
%               Triangle, Green: Time-Position of the HRIR-Maximum
%               Diamonds, Blue: Centroid of the HRIR
%               Squares, Magenta: Average Group Delay (1-5 kHz)
%               
%               Right panel: 
%               Estimated TOAs of NH89 (ARI database) in the horizontal
%               interaural plane.
%               Black: Minimum-Phase Cross-Correlation Method
%               Blue: Time-Position of the HRIR-Maximum
%               Green: Centroid of the HRIR
%               Red: Average Group Delay (1-5 kHz)
%
%     'fig3'    Reproduce Fig. 3:
%               
%               Left panel: 
%               Sagittal TOA deviations and averaged TOA variance for NH89
%               as function of the polar angle for all sagittal groups.
%               Dots, Blue: Sagittal TOA deviations
%               Line, Red: Averaged TOA variance
%               
%               Right panel: 
%               Estimated TOAs, detected outliers and outlier adjusted set
%               of TOAs for NH89 (ARI database) in the horizontal plane.
%               Line: Estimated TOAs
%               Triangles (down), Blue: Detected outliers for the azimuthal
%               slope criterion
%               Triangles (up), Blue: Detected outliers for the sagittal
%               variance criterion
%               Dots, Red: Outlier-adjusted set
%
%     'fig5'    Reproduce Fig. 5:
%               Model parameters (sphere radius, ear position) resulting
%               from fitting the on-axis model to HRTFs of a rigid sphere.
%               Squares: for the left ear
%               Diamonds: for the right ear
%               Dashed lines: set values used in the numerical HRTF
%               calculation
%
%     'fig6'    Reproduce Fig. 6:
%               TOA in interaural horizontal plane for the left ear HRTFs
%               of NH89.
%               Solid line, Black: On-axis model fitted to outlier-adjusted
%               set of TOAs
%               Circles, Red: Outlier-adjusted set of TOAs
%               Hexagrams, Blue: Detected outliers
%
%     'fig7'    Reproduce Fig. 7:
%               Model parameter (sphere radius) resulting from fitting the
%               on-axis model to HRTFs of human listeners. The listeners
%               are sorted by the ascending binaural average radius.
%               Blue: for the left ear
%               Green: for the right ear
%               Circles: ARI
%               Diamonds: CIPIC
%               Squares: LISTEN
%
%     'fig8'    Reproduce Fig. 8:
%               Relative TOAs for NH89 (ARI database) in the interaural
%               horizontal plane.
%               Dashed lines: for the right ear
%               Solid lines: for the right ear
%               Thin lines: TOAs estimated with the Minimum-phase
%               Cross-Correlation method
%               Thick lines: On-axis model fitted to the outlier-adjusted
%               sets of TOAs
%
%     'fig9'    Reproduce Fig. 9:
%               Model parameters (sphere radius) resulting from fitting the
%               on-axis model to HRTFs of an off-axis placed rigid sphere.
%               All other conventions as in Fig. 5.
%
%     'fig11'   Reproduce Fig. 11:
%               Model parameters (sphere radius, sphere center) resulting
%               from fitting the off-axis model to HRTFs of a rigid
%               sphere.
%               All other conventions as in Fig. 5.
%
%     'fig12'   Reproduce Fig. 12:
%               Model parameters (sphere radius, sphere center) resulting
%               from fitting the off-axis model to HRTFs of human
%               listeners.
%               All other conventions as in Fig. 7.
%
%   Examples:
%   ---------
%
%   To display Fig. 2, use :::
%
%     exp_ziegelwanger2014('fig2');
%
%   To display Fig. 3, use :::
%
%     exp_ziegelwanger2014('fig3');
%
%   To display Fig. 5, use :::
%
%     exp_ziegelwanger2014('fig5');
%
%   To display Fig. 7, use :::
%
%     exp_ziegelwanger2014('fig7');
%
%   To display Fig. 8, use :::
%
%     exp_ziegelwanger2014('fig8');
%
%   To display Fig. 9, use :::
%
%     exp_ziegelwanger2014('fig9');
%
%   To display Fig. 11, use :::
%
%     exp_ziegelwanger2014('fig11');
%
%   To display Fig. 12, use :::
%
%     exp_ziegelwanger2014('fig12');
%
%   See also: ziegelwanger2014, ziegelwanger2014onaxis,
%   ziegelwanger2014offaxis, data_ziegelwanger2014
%
%   References: ziegelwanger2014

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria

%% ------ Check input options --------------------------------------------

    definput.flags.type = {'missingflag',...
    'fig2','fig3','fig5','fig6',...
    'fig7','fig8','fig9','fig11','fig12','fig5new','fig6new','fig7new','fig8new','fig10new','fig12new','tab2'};
    definput.flags.plot = {'plot','noplot'};
    definput.flags.results = {'reload','recalc'};

    % Parse input options
    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
        flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}), ...
            sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
        error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
    end;

%% Figure 2
if flags.do_fig2
    
    data=data_ziegelwanger2014('NH89');

    subplot(122)
    %---------------------------Threshold---------------------------
    [~,tmp]=ziegelwanger2014(data,1,0,0);
    MAX=tmp.toa;

    %---------------------------Centroid----------------------------
    [~,tmp]=ziegelwanger2014(data,2,0,0);
    CTD=tmp.toa;

    %---------------------------Groupdelay--------------------------
    [~,tmp]=ziegelwanger2014(data,3,0,0);
    AGD=tmp.toa;

    %---------------------------Minimal-Phase-----------------------
    [~,tmp]=ziegelwanger2014(data,4,0,0);
    MCM=tmp.toa;
    clear tmp

    plotziegelwanger2014(data,MCM(:,1),4,[0 0 0]/255,0,1,1,'-',1);
    hold on
    plotziegelwanger2014(data,MAX(:,1),4,[0 0 80]/255,0,1,1,'--',1);
    plotziegelwanger2014(data,CTD(:,1),4,[50 220 50]/255,0,1,1,'-',1);
    plotziegelwanger2014(data,AGD(:,1),4,[250 80 80]/255,0,1,1,'--',1);
    xlim([-10 370])
    ylim([2.65 4.05])
    grid off
    legend(' MCM',' MAX',' CTD',' AGD','Location','NorthWest');
    legend boxoff
    xlabel('Azimuth (in deg) ')
    ylabel('TOA (ms) ')
%     set(gca,'FontSize',ticklabelsize,'LineWidth',bw,'LineWidth',bw);%,'YTick',[-0.8 -0.4 0 0.4 0.8]
    title('');
    
    subplot(121)
    time=(0:data.DimSize.N-1)/data.Data.SamplingRate*1000;
    MAX=round(MAX);
    CTD=round(CTD);
    AGD=round(AGD);
    MCM=round(MCM);
    idx=ARI_FindPosition(data,85,0);

    fprintf(['MAX: ITD is ' num2str(time(diff(MAX(idx,:),1,2))) ' ms\n'])
    fprintf(['CTD: ITD is ' num2str(time(diff(CTD(idx,:),1,2))) ' ms\n'])
    fprintf(['AGD: ITD is ' num2str(time(diff(AGD(idx,:),1,2))) ' ms\n'])
    fprintf(['MCM: ITD is ' num2str(time(diff(MCM(idx,:),1,2))) ' ms\n'])

    plot(time,squeeze(data.Data.IR(idx,1,:))/max(abs(data.Data.IR(idx,1,:))),'k-');
    hold on
    plot(time,squeeze(data.Data.IR(idx,2,:))/max(abs(data.Data.IR(idx,2,:))),'k--')
    h=stem([time(MCM(idx,1)) time(MCM(idx,1))],[-2 data.Data.IR(idx,1,MCM(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'r-','BaseValue',-1);
    stem([time(CTD(idx,1)) time(CTD(idx,1))],[-2 data.Data.IR(idx,1,CTD(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'g-','BaseValue',-1)
    stem([time(MAX(idx,1)) time(MAX(idx,1))],[-2 data.Data.IR(idx,1,MAX(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'b-','BaseValue',-1)
    stem([time(AGD(idx,1)) time(AGD(idx,1))],[-2 data.Data.IR(idx,1,AGD(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'m-','BaseValue',-1)
    stem([time(MCM(idx,2)) time(MCM(idx,2))],[-2 data.Data.IR(idx,2,MCM(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'r-','BaseValue',-1)
    stem([time(CTD(idx,2)) time(CTD(idx,2))],[-2 data.Data.IR(idx,2,CTD(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'g-','BaseValue',-1)
    stem([time(AGD(idx,2)) time(AGD(idx,2))],[-2 data.Data.IR(idx,2,AGD(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'m-','BaseValue',-1)
    stem([time(MAX(idx,2)) time(MAX(idx,2))],[-2 data.Data.IR(idx,2,MAX(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'b-','BaseValue',-1)
    plot(time(MAX(idx,1)),data.Data.IR(idx,1,MAX(idx,1))/max(abs(data.Data.IR(idx,1,:))),'b^','MarkerFaceColor','b')
    plot(time(MCM(idx,1)),data.Data.IR(idx,1,MCM(idx,1))/max(abs(data.Data.IR(idx,1,:))),'ro','MarkerFaceColor','r')
    plot(time(CTD(idx,1)),data.Data.IR(idx,1,CTD(idx,1))/max(abs(data.Data.IR(idx,1,:))),'gd','MarkerFaceColor','g')
    plot(time(AGD(idx,1)),data.Data.IR(idx,1,AGD(idx,1))/max(abs(data.Data.IR(idx,1,:))),'ms','MarkerFaceColor','m')
    plot(time(MAX(idx,2)),data.Data.IR(idx,2,MAX(idx,2))/max(abs(data.Data.IR(idx,2,:))),'b^','MarkerFaceColor','b')
    plot(time(MCM(idx,2)),data.Data.IR(idx,2,MCM(idx,2))/max(abs(data.Data.IR(idx,2,:))),'ro','MarkerFaceColor','r')
    plot(time(CTD(idx,2)),data.Data.IR(idx,2,CTD(idx,2))/max(abs(data.Data.IR(idx,2,:))),'gd','MarkerFaceColor','g')
    plot(time(AGD(idx,2)),data.Data.IR(idx,2,AGD(idx,2))/max(abs(data.Data.IR(idx,2,:))),'ms','MarkerFaceColor','m')
    xlim([2.4 4.1])
    ylim([-1.1 1.1])
    xlabel('Time (ms) ')
    ylabel('Amplitude ')
%     set(gca,'FontSize',ticklabelsize,'LineWidth',bw,'LineWidth',bw,'XTick',[2.5 3 3.5 4]);
    set(get(h,'Baseline'),'Visible','off')
end

%% Figure 3, Figure 6
if flags.do_fig3 || flags.do_fig6
        
    data=data_ziegelwanger2014('NH89');

    p0_onaxis=[[0.0875; pi/2; 0; 0.0001] [0.0875; -pi/2; 0; 0.0001]];
    p0_onaxis=transpose(p0_onaxis);
    p_onaxis=zeros(size(p0_onaxis));
    p0_offaxis=zeros(2,7);
    p_offaxis=p0_offaxis;

    toa=zeros(data.DimSize.M,data.DimSize.R);
    toaEst=zeros(data.DimSize.M,data.DimSize.R);
    indicator=zeros(data.DimSize.M,data.DimSize.R);
    indicator_hor=indicator;
    indicator_sag=indicator;
    pos=zeros(data.DimSize.M,8);
    pos(:,1:2)=data.APV(:,1:2);
    [pos(:,6),pos(:,7)]=sph2hor(data.APV(:,1),data.APV(:,2));
    pos(:,8)=cumsum(ones(data.DimSize.M,1));
    [~,tmp]=ziegelwanger2014(data,4,0,0);
    toaEst=tmp.toa;
    
    for ch=1:data.DimSize.R

        % Outlier detection: smooth TOA in horizontal planes
        [~,idxSortHor]=sort(pos(:,1));
        epsilon=5;
        slope=zeros(data.DimSize.M,1);
        for ele=min(pos(:,2)):epsilon:max(pos(:,2)) %calculate slope for each elevation along azimuth
            idx=find(pos(idxSortHor,2)>ele-epsilon/2 & pos(idxSortHor,2)<=ele+epsilon/2);
            if numel(idx)>1
                idx(length(idx)+1)=idx(1);
                slope(idxSortHor(idx(1:end-1)),1)=diff(toaEst(idxSortHor(idx),ch))./(abs(abs(abs(diff(pos(idxSortHor(idx),1)))-180)-180)+0.00000000001);
            end
        end
        sloperms=sqrt(sum(slope.^2)/length(slope));
        if sloperms<30/(length(find(pos(:,2)==0))/2)
            sloperms=30/(length(find(pos(:,2)==0))/2);
        end
        for ele=min(pos(:,2)):epsilon:max(pos(:,2))
            idx=find(pos(idxSortHor,2)>ele-epsilon/2 & pos(idxSortHor,2)<=ele+epsilon/2);
            for ii=1:length(idx)-1
                if abs(slope(idxSortHor(idx(ii))))>sloperms
                    for jj=0:1
                        if ii+jj==length(idx)
                            indicator_hor(idxSortHor(idx(end)),ch)=1;
                        else
                            indicator_hor(idxSortHor(idx(mod(ii+jj,length(idx)))),ch)=1;
                        end
                    end
                end
            end
            clear idx
        end

        % Outlier detection: constant TOA in sagittal planes
        epsilon=2;
        for ii=1:20
            sag_dev=zeros(data.DimSize.M,1);
            for lat=-90:epsilon:90
                idx=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2); 
                idx2=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2 & indicator_hor(:,ch)==0 & indicator(:,ch)==0);
                if length(idx2)>2
                    sag_dev(idx,1)=toaEst(idx,ch)-mean(toaEst(idx2,ch));
                end
            end
            sag_var=sqrt(sum(sag_dev.^2)/length(sag_dev));
            if sag_var<2
                sag_var=2;
            end
            indicator(:,ch)=zeros(size(indicator,1),1);
            indicator_sag(:,ch)=zeros(size(indicator_sag,1),1);
            indicator_sag(abs(sag_dev)>sag_var,ch)=ones(length(find(abs(sag_dev)>sag_var)),1);
            indicator(abs(sag_dev)>sag_var | indicator_hor(:,ch)==1,ch)=ones(length(find(abs(sag_dev)>sag_var | indicator_hor(:,ch)==1)),1);
        end

        if flags.do_fig3 && ch==1 %Figure 3
            subplot(121)
            plot([-90; 270],[sag_var/data.Data.SamplingRate*1000000; sag_var/data.Data.SamplingRate*1000000],'r--');
            hold on
            plot(real(pos(:,7)),abs(sag_dev)/data.Data.SamplingRate*1000000,'b.');
            xlim([-98 278])
            ylim([-5 max(abs(sag_dev)/data.Data.SamplingRate*1000000)+5])
            xlabel('Polar angle in degree')
            ylabel('Sagittal TOA deviation in Âµs')
            title('')
            set(gca,'XTick',[-90 0 90 180 270])
        end
        clear sag_dev; clear sag_var;

        if flags.do_fig3 && ch==1 %Figure 3
            subplot(122)
            plotziegelwanger2014(data,toaEst,4,'k',0,1,1,{'-'},1);
            hold on
            h=plotziegelwanger2014(data,indicator_sag.*toaEst,3,'w',0,1,1,{'^'},4);
            set(h,'MarkerFaceColor','w','MarkerEdgeColor','w');
            h=plotziegelwanger2014(data,indicator_hor.*toaEst,3,'w',0,1,1,{'v'},4);
            set(h,'MarkerFaceColor','w','MarkerEdgeColor','w');
            h=plotziegelwanger2014(data,indicator_sag.*toaEst,3,'b',0,1,1,{'^'},4);
            set(h,'LineWidth',2);
            h=plotziegelwanger2014(data,indicator_hor.*toaEst,3,'b',0,1,1,{'v'},4);
            set(h,'LineWidth',2);
            h=plotziegelwanger2014(data,(-indicator+1).*toaEst,3,'r',0,1,1,{'o'},4);
            set(h,'MarkerFaceColor','r','MarkerEdgeColor','r');
            ylabel('TOA in ms')
            xlabel('Azimuth in degree')
            grid off
            xlim([-10 370])
            ylim([2.65 3.65])
            title('')
            set(gca,'YTick',[2.7 3 3.3 3.6])
        end
    end

    for ch=1:data.DimSize.R
        p0_onaxis(ch,4)=min(toaEst(indicator(:,ch)==0,ch))/data.Data.SamplingRate;
        p0offset_onaxis=[0.06 pi/4 pi/4 0.001];

        idx=find(indicator(:,ch)==0);
        x=pos(idx,1:2)*pi/180;
        y=toaEst(idx,ch)/data.Data.SamplingRate;
        p_onaxis(ch,:)=lsqcurvefit(@ziegelwanger2014onaxis,p0_onaxis(ch,:),x,y,p0_onaxis(ch,:)-p0offset_onaxis,p0_onaxis(ch,:)+p0offset_onaxis,optimset('Display','off','TolFun',1e-6));
        toa(:,ch)=ziegelwanger2014onaxis(p_onaxis(ch,:),pos(:,1:2)*pi/180)*data.Data.SamplingRate;
    end

    TolFun=[1e-5; 1e-6];
    for ii=1:size(TolFun,1)
        for ch=1:data.DimSize.R
            idx=find(indicator(:,ch)==0);
            x=pos(idx,1:2)*pi/180;
            y=toaEst(idx,ch)/data.Data.SamplingRate;
            p0_offaxis(ch,:)=[p0_onaxis(ch,1) 0 0 0 p0_onaxis(ch,4) p0_onaxis(ch,2) p0_onaxis(ch,3)];
            p0offset_offaxis=[0.05 0.05 0.05 0.05 0.001 pi pi];
            p_offaxis(ch,:)=lsqcurvefit(@ziegelwanger2014offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0offset_offaxis,p0_offaxis(ch,:)+p0offset_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));
            toa(:,ch)=ziegelwanger2014offaxis(p_offaxis(ch,:),pos(:,1:2)*pi/180)*data.Data.SamplingRate;
        end
        if abs(diff(p_offaxis(:,1)))>0.003 || abs(diff(p_offaxis(:,3)))>0.003
            p_offaxis(:,[1 3])=p_offaxis([2 1],[1 3]);
            for ch=1:data.DimSize.R
                idx=find(indicator(:,ch)==0);
                x=pos(idx,1:2)*pi/180;
                y=toaEst(idx,ch)/data.Data.SamplingRate;
                p0_offaxis(ch,:)=[p_offaxis(ch,1) mean(p_offaxis(:,2)) p_offaxis(ch,3) mean(p_offaxis(:,4)) mean(p_offaxis(:,5)) p_offaxis(ch,6) p_offaxis(ch,7)];
                p0offset_offaxis=[0.05 0.05 0.05 0.05 0.001 pi/2 pi/2];
                p_offaxis(ch,:)=lsqcurvefit(@ziegelwanger2014offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0offset_offaxis,p0_offaxis(ch,:)+p0offset_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));
                toa(:,ch)=ziegelwanger2014offaxis(p_offaxis(ch,:),pos(:,1:2)*pi/180)*data.Data.SamplingRate;
            end
        end
        if abs(diff(p_offaxis(:,1)))<0.003 && abs(diff(p_offaxis(:,2)))<0.003 && abs(diff(p_offaxis(:,3)))<0.003 && abs(diff(p_offaxis(:,4)))<0.003
            break
        end
    end

    if flags.do_fig6 %Figure 6
        h=plotziegelwanger2014(data,indicator.*toaEst,3,'b',0,1,1,{'^'},4);
        set(h,'LineWidth',2);
        hold on
        h=plotziegelwanger2014(data,indicator.*toaEst,3,'b',0,1,1,{'v'},4);
        set(h,'LineWidth',2);
        h=plotziegelwanger2014(data,(-indicator+1).*toaEst,3,'r',0,1,1,{'o'},4);
        set(h,'MarkerFaceColor','r','MarkerEdgeColor','r');
        plotziegelwanger2014(data,toa,4,'k',0,1,1,{'-'},1);
        ylabel('TOA ms ')
        xlabel('Azimuth in degree')
        grid off
        xlim([-10 370])
        ylim([2.65 3.65])
        title('')
        set(gca,'YTick',[2.7 3 3.3 3.6])
    end
end

%% Figure 5
if flags.do_fig5

    sym='sdo'; %plot symbols
    clr=[0,0,255; 255,0,0; 255,255,67]/255; %plot colors
    meclr=[0,0,255; 255,0,0; 255,255,67]/255; %marker edge colors
    
    if flags.do_recalc
        data=data_ziegelwanger2014('SPHERE_ROT','recalc');
    else
        data=data_ziegelwanger2014('SPHERE_ROT');
    end
    
    % radii
    subplot(311)
    var=[squeeze(data.results.p_onaxis(1,1,:))*100 squeeze(data.results.p_onaxis(1,2,:))*100 data.radius(:)/10];
    err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
    err=reshape(err,numel(err),1);
    fprintf(['Radius: average err is ' num2str(mean(err)) ' cm\n']);
    fprintf(['        standard deviation is ' num2str(std(err)) ' cm\n']);
    for ch=1:size(data.results.p_onaxis,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    plot(var(:,3),'k--')
    clear var;
    ylabel('r in cm')

    %phi
    subplot(312)
    var=[squeeze(data.results.p_onaxis(2,1,:))/pi*180 squeeze(data.results.p_onaxis(2,2,:))/pi*180 -data.phi+ones(length(data.phi),1)*90 -data.phi-ones(length(data.phi),1)*90];
    err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,4)]);
    err=reshape(err,numel(err),1);
    fprintf(['Phi: average err is ' num2str(mean(err)) 'deg\n']);
    fprintf(['     standard deviation is ' num2str(std(err)) 'deg\n']);
    for ch=1:size(data.results.p_onaxis,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    for ch=1:size(data.results.p_onaxis,2)
        plot(1:9,var(1:9,2+ch),'k--')
        plot(10:14,var(10:14,2+ch),'k--')
        plot(15:23,var(15:23,2+ch),'k--')
        plot(24:28,var(24:28,2+ch),'k--')
        plot(29:37,var(29:37,2+ch),'k--')
        plot(38:42,var(38:42,2+ch),'k--')
    end
    clear var;
    set(gca,'YTick',[-90 90])
    ylabel('\phi_e in deg')

    %theta
    subplot(313)
    var=[squeeze(data.results.p_onaxis(3,1,:))/pi*180 squeeze(data.results.p_onaxis(3,2,:))/pi*180 data.theta -data.theta];
    err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,4)]);
    err=reshape(err,numel(err),1);
    fprintf(['Theta: average err is ' num2str(mean(err)) 'deg\n']);
    fprintf(['       standard deviation is ' num2str(std(err)) 'deg\n']);
    for ch=1:size(data.results.p_onaxis,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    for ch=1:size(data.results.p_onaxis,2)
        plot(1:9,var(1:9,2+ch),'k--')
        plot(10:14,var(10:14,2+ch),'k--')
        plot(15:23,var(15:23,2+ch),'k--')
        plot(24:28,var(24:28,2+ch),'k--')
        plot(29:37,var(29:37,2+ch),'k--')
        plot(38:42,var(38:42,2+ch),'k--')
    end
    clear var;
    ylabel('\theta_e in deg')
    xlabel('Condition')
end

%% Figure 7, Figure 12
if flags.do_fig7 || flags.do_fig12
    
    hrtf={'ARI','CIPIC','LISTEN'};
    sym='ods'; %plot symbols

    %-------------------------------Load Data----------------------------------
    for kk=1:length(hrtf)
        if flags.do_recalc
            data=data_ziegelwanger2014(hrtf{kk},'recalc');
        else
            data=data_ziegelwanger2014(hrtf{kk});
        end
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        temp1=zeros(size(data.results(1).meta.p_onaxis,1),size(data.results(1).meta.p_onaxis,2),length(data.results));
        temp3=zeros(size(data.results(1).meta.p_offaxis,1),size(data.results(1).meta.p_offaxis,2),length(data.results));
        temp4=zeros(length(data.results),4);
        temp6=zeros(length(data.results),4);
        temp7=zeros(length(data.results),4);
        for ii=1:length(data.results)
            temp1(:,:,ii)=data.results(ii).meta.p_onaxis;
            temp3(:,1:size(data.results(ii).meta.p_offaxis,2),ii)=data.results(ii).meta.p_offaxis;
            temp4(ii,:)=[data.results(ii).meta.performanceOutliers(1).outlierRate data.results(ii).meta.performanceOutliers(2).outlierRate data.results(ii).meta.performanceOutliers(3).outlierRate data.results(ii).meta.performanceOutliers(4).outlierRate];
            if isfield(data.results(ii).meta.performanceOutliers(1)','outliersl')
                temp6(ii,:)=[data.results(ii).meta.performanceOutliers(1).outlierRateL data.results(ii).meta.performanceOutliers(2).outlierRateL data.results(ii).meta.performanceOutliers(3).outlierRateL data.results(ii).meta.performanceOutliers(4).outlierRateL];
                temp7(ii,:)=[data.results(ii).meta.performanceOutliers(1).outlierRateR data.results(ii).meta.performanceOutliers(2).outlierRateR data.results(ii).meta.performanceOutliers(3).outlierRateR data.results(ii).meta.performanceOutliers(4).outlierRateR];
            end
        end
        p_onaxis{kk}=temp1;
        p_offaxis{kk}=temp3;
        outlierRate{kk}=temp4;
        outliersear{kk}=[temp6; temp7];
    end
    
    if flags.do_fig7 %Figure 7
        fprintf('On-Axis Model:\n')
        fprintf(['Average Radius: ' ...
                 num2str(mean([mean(squeeze(p_onaxis{1}(1,:,:)*100)) ...
                            mean(squeeze(p_onaxis{2}(1,:,:)*100)) ...
                            mean(squeeze(p_onaxis{3}(1,:,:)*100))])) ' cm\n'])
        fprintf(['Average Radius Difference: ' num2str(mean([mean(abs(diff(squeeze(p_onaxis{1}(1,:,:)*100)))) mean(abs(diff(squeeze(p_onaxis{2}(1,:,:)*100)))) mean(abs(diff(squeeze(p_onaxis{3}(1,:,:)*100))))])) 'cm\n'])
        fprintf(['Maximum Radius Difference: ' num2str(max([max(abs(diff(squeeze(p_onaxis{1}(1,:,:)*100)))) max(abs(diff(squeeze(p_onaxis{2}(1,:,:)*100)))) max(abs(diff(squeeze(p_onaxis{3}(1,:,:)*100))))])) ' cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(mean(p_onaxis{kk}(1,:,:)*100,2)) kk* ...
                 ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                        size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(p_onaxis{kk}(1,1,:)*100) kk*ones(size(p_onaxis{kk},3),1) ...
                 transpose(1:size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(p_onaxis{kk}(1,2,:)*100) kk*ones(size(p_onaxis{kk},3),1) ...
                 transpose(1:size(p_onaxis{kk},3))];
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
            h{kk}=plot(varl(varl(:,2)==kk,3),varl(varl(:,2)==kk,1), ...
                       sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
            plot(varr(varr(:,2)==kk,3),varr(varr(:,2)==kk,1),sym(kk), ...
                 'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        clear var varl varr stm;
        xlabel('Listeners ')
        ylabel('r in cm')
        xlim([-1.5 temp+1.5])
        ylim([4.1 12+0.7])
        legend([h{1},h{2},h{3}],'ARI','CIPIC','LISTEN','Location','NorthWest');
        legend boxoff
    end
    
    if flags.do_fig12 %Figure12
        %radii
        subplot(411)
        fprintf('Off-Axis Model:\n')
        fprintf(['Radius (Avg): ' ...
                 num2str(mean([mean(squeeze(p_offaxis{1}(1,:,:)*100)) ...
                            mean(squeeze(p_offaxis{2}(1,:,:)*100)) ...
                            mean(squeeze(p_offaxis{3}(1,:,:)*100))])) ' cm\n'])
        fprintf(['Radius (Std): ' num2str(std([reshape(p_offaxis{1}(1,:,: ...
                                                          )*100,1, ...
                                                       numel(p_offaxis{1}(1,:,:))) reshape(p_offaxis{2}(1,:,:)*100,1,numel(p_offaxis{2}(1,:,:))) reshape(p_offaxis{3}(1,:,:)*100,1,numel(p_offaxis{3}(1,:,:)))])) ' cm\n'])
        fprintf(['Radius Difference (Avg): ' ...
                 num2str(mean([mean(abs(diff(squeeze(p_offaxis{1}(1,:,:)*100)))) ...
                            mean(abs(diff(squeeze(p_offaxis{2}(1,:,:)*100)))) ...
                            mean(abs(diff(squeeze(p_offaxis{3}(1,:,:)*100))))])) ...
                 'cm\n'])
        fprintf(['Radius Difference (Std): ' ...
                 num2str(std([abs(diff(squeeze(p_offaxis{1}(1,:,:)*100))) ...
                            abs(diff(squeeze(p_offaxis{2}(1,:,:)*100))) ...
                            abs(diff(squeeze(p_offaxis{3}(1,:,:)*100)))])) ...
                 'cm\n'])
        fprintf(['Radius Difference (Max): ' ...
                 num2str(max([max(abs(diff(squeeze(p_offaxis{1}(1,:,:)*100)))) ...
                            max(abs(diff(squeeze(p_offaxis{2}(1,:,:)*100)))) ...
                            max(abs(diff(squeeze(p_offaxis{3}(1,:,:)*100))))])) ...
                 ' cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(mean(p_offaxis{kk}(1,:,:)*100,2)) kk* ...
                 ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                        size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(p_offaxis{kk}(1,1,:)*100) kk* ...
                 ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                        size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(p_offaxis{kk}(1,2,:)*100) kk* ...
                 ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                        size(p_onaxis{kk},3))];
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
        ylabel('r in cm')
        xlim([-1.5 temp+1.5])
        ylim([4.1 12+0.7])
        legend([h{1},h{2},h{3}],'ARI','CIPIC','LISTEN','Location','NorthWest');
        legend boxoff

        %lateral displacements
        % xM
        subplot(412)
        fprintf(['Average xM: ' ...
                 num2str(mean([mean(abs(squeeze(p_offaxis{1}(2,:,:)*100))) ...
                            mean(abs(squeeze(p_offaxis{2}(2,:,:)*100))) ...
                            mean(abs(squeeze(p_offaxis{3}(2,:,:)*100)))])) ...
                 ' cm\n'])
        fprintf(['Average xM Difference: ' ...
                 num2str(mean([mean(abs(diff(squeeze(p_offaxis{1}(2,:,:)*100)))) ...
                            mean(abs(diff(squeeze(p_offaxis{2}(2,:,:)*100)))) ...
                            mean(abs(diff(squeeze(p_offaxis{3}(2,:,:)*100))))])) ...
                 'cm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(mean(p_offaxis{kk}(2,:,:)*100,2)) kk* ...
                 ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                        size(p_onaxis{kk},3))];
            varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(p_offaxis{kk}(2,1,:)*100) kk* ...
                 ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                        size(p_onaxis{kk},3))];
            varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
                [squeeze(p_offaxis{kk}(2,2,:)*100) kk* ...
                 ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                        size(p_onaxis{kk},3))];
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
        ylabel('x_M in cm')
        xlim([-1.5 temp+1.5])
        ylim([-4.5 4.5])

        % yM
        subplot(413)
        fprintf(['yM (Avg): ' num2str(mean([mean(abs(squeeze(p_offaxis{1}(3,:,:)*100))) mean(abs(squeeze(p_offaxis{2}(3,:,:)*100))) mean(abs(squeeze(p_offaxis{3}(3,:,:)*100)))])) ' cm\n'])
        fprintf(['yM (Std): ' num2str(std([reshape(p_offaxis{1}(3,:,:)*100,1,numel(p_offaxis{1}(3,:,:))) reshape(p_offaxis{2}(3,:,:)*100,1,numel(p_offaxis{2}(3,:,:))) reshape(p_offaxis{3}(3,:,:)*100,1,numel(p_offaxis{3}(3,:,:)))])) ' cm\n'])
        fprintf(['yM Difference (Avg): ' num2str(mean([mean(abs(diff(squeeze(p_offaxis{1}(3,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{2}(3,:,:)*100)))) mean(abs(diff(squeeze(p_offaxis{3}(3,:,:)*100))))])) 'cm\n'])
        fprintf(['yM Difference (Std): ' num2str(std([abs(diff(squeeze(p_offaxis{1}(3,:,:)*100))) abs(diff(squeeze(p_offaxis{2}(3,:,:)*100))) abs(diff(squeeze(p_offaxis{3}(3,:,:)*100)))])) 'cm\n'])
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
        ylabel('y_M in cm')
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
        ylabel('z_M in cm')
        xlim([-1.5 temp+1.5])
        ylim([-4.5 6.5])
    end
    
end

%% Figure 8
if flags.do_fig8
    
    Obj=data_ziegelwanger2014('NH89');
    
    [~,tmp]=ziegelwanger2014(Obj,4,1);
    toa1=tmp.toa;
    fprintf(['Radii for left and right ear: ' num2str(tmp.p_onaxis(1,:)*100) ' cm\n'])
    [~,tmp]=ziegelwanger2014(Obj,4,0,0);
    toa2=tmp.toa;
    fprintf(['Maximum TOA difference left: ' num2str((max(toa1(:,1))-min(toa1(:,1)))/Obj.Data.SamplingRate*1000) ' ms\n'])
    fprintf(['Maximum TOA difference right: ' num2str((max(toa2(:,2))-min(toa1(:,2)))/Obj.Data.SamplingRate*1000) ' ms\n'])

    plotziegelwanger2014(Obj,toa2-min(min(toa1)),4,'k',0,1,1,'--',1);
    plotziegelwanger2014(Obj,toa1-min(min(toa1)),4,'k',0,1,1,'--',2);
    plotziegelwanger2014(Obj,toa2-min(min(toa1)),4,'k',0,2,1,'-',1);
    plotziegelwanger2014(Obj,toa1-min(min(toa1)),4,'k',0,2,1,'-',2);

    xlim([-5 365])
    ylim([-0.05 0.95])
    grid off
    xlabel('Azimuth in degree')
    ylabel('Relative TOA in ms');
    title('')
    
end

%% Figure 9, Figure 11
if flags.do_fig9 || flags.do_fig11
    
    sym='sdo'; %plot symbols
    clr=[0,0,255; 255,0,0; 255,255,67]/255; %plot colors
    meclr=[0,0,255; 255,0,0; 255,255,67]/255; %marker edge colors
    
    if flags.do_recalc
        data=data_ziegelwanger2014('SPHERE_DIS','recalc');
    else
        data=data_ziegelwanger2014('SPHERE_DIS');
    end
    
    if flags.do_fig9 %Figure 9
        p1=data.results.p_onaxis(:,:,[1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);
        r1=data.radius([1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);
        yM1=data.yM([1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);
        
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
        ylabel('r in cm')

        %yM
        subplot(212)
        plot(-yM1*100,'k--')
        clear var;
        xlabel('Condition')
        ylabel('y_M in cm ')
    end
    
    if flags.do_fig11 %Figure 11
        center=[data.xM(1:length(data.xM)/3) data.yM(1:length(data.yM)/3) data.zM(1:length(data.zM)/3)];
        [~,idx3]=sort(squeeze(center(:,3)));
        [~,idx2]=sort(squeeze(center(:,1)));
        [~,idx1]=sort(squeeze(center(:,2)));
        idx=idx3(idx2(idx1));
        idx=[idx; idx+length(data.xM)/3; idx+length(data.xM)/3*2];
        data.radius=data.radius(idx);
        p_offaxis=data.results.p_offaxis(:,:,idx);
        data.xM=data.xM(idx);
        data.yM=data.yM(idx);
        data.zM=data.zM(idx);

        %radii
        subplot(411)
        ymax2=round(max(max(squeeze(p_offaxis(1,:,:)*100))))+1;
        var=[squeeze(p_offaxis(1,1,:))*100 squeeze(p_offaxis(1,2,:))*100 data.radius/10];
        err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
        err=reshape(err,numel(err),1);
        temp1=err;
        fprintf(['Radius: average err is ' num2str(mean(err)) ' cm\n']);
        fprintf(['        standard deviation is ' num2str(std(err)) ' cm\n']);
        fprintf(['        maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
        fprintf(['        average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        for ch=1:size(p_offaxis,2)
            plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
            hold on
        end
        plot(var(:,3),'k--')
        ylabel('r in cm')
        xlim([0 size(var,1)+1])
        ylim([4.5 ymax2])
        clear var;

        %xM
        subplot(412)
        var=[squeeze(p_offaxis(2,1,:))*100 squeeze(p_offaxis(2,2,:))*100 -data.xM*100];
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
        ylabel('x_M in cm ')
        xlabel([0 size(var,1)+1])
        xlim([0 size(var,1)+1])
        ylim([-2.5 0.5])
        clear var;

        %yM
        subplot(413)
        var=[squeeze(p_offaxis(3,1,:))*100 squeeze(p_offaxis(3,2,:))*100 -data.yM*100];
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
        ylabel('y_M in cm ')
        xlim([0 size(var,1)+1])
        ylim([-0.5 2.5])
        clear var;

        %zM
        subplot(414)
        var=[squeeze(p_offaxis(4,1,:))*100 squeeze(p_offaxis(4,2,:))*100 -data.zM*100];
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
        set(gca,'YTick',[-1 0])
        xlabel('Condition')
        ylabel('z_M in cm')
        xlim([0 size(var,1)+1])
        ylim([-1.5 0.5])
        clear var;

        
        fprintf(['offset: average err is ' num2str(mean(temp1)) ' cm\n']);
        fprintf(['        standard deviation is ' num2str(std(temp1)) ' cm\n']);
        fprintf(['        maximum is ' num2str(max(abs(temp1))) ' cm\n']);
    end
end

%% Figure 5 new
if flags.do_fig5new
    cd '/Volumes/ARI/Simulations/HRTF-TOA/SphereTorsoPinna/'
    method=1;
    methodLabel=['MAX';'CTD';'AGD';'MCM'];
    clr='kbgr';
    
    load(['Sphere' filesep 'hrtf_M_hrtf_2ears.mat']);
    Obj1=SOFAconvertARI2SOFA(hM,meta,stimPar);
    load(['SphereTorso' filesep 'hrtf_M_hrtf_2ears.mat']);
    Obj2=SOFAconvertARI2SOFA(hM,meta,stimPar);
    load(['SphereTorsoPinna' filesep 'hrtf_M_hrtf_2ears.mat']);
    Obj3=SOFAconvertARI2SOFA(hM,meta,stimPar);

    figure
    [Obj1,~]=ziegelwanger2014(Obj1,4,0,0);
    plotziegelwanger2014(Obj1,Obj1.Data.Delay,1,clr(1),0,1,0);
    hold on
    [Obj2,~]=ziegelwanger2014(Obj2,4,0,0);
    plotziegelwanger2014(Obj2,Obj2.Data.Delay,1,clr(2),0,1,0);
    [Obj3,~]=ziegelwanger2014(Obj3,4,0,0);
    plotziegelwanger2014(Obj3,Obj3.Data.Delay,1,clr(3),0,1,0);
    legend('SPH','SPH+TOR','SPH+TOR+PIN','Location','NorthWest')
    export_fig('~/Dropbox/HRTF-TOA (1)/Publications/Paper (Model)/2 Revision 1/Figures New/Fig5l.png','-png','-transparent','-m2')
    close

    figure
    for method=1:4
        [Obj3,results]=ziegelwanger2014(Obj3,method,0,0);
        plotziegelwanger2014(Obj3,Obj3.Data.Delay,1,clr(method),0,1,0);
        hold on
    end
    legend('MAX','CTD','AGD','MCM','Location','NorthWest')
    
    export_fig('~/Dropbox/HRTF-TOA (1)/Publications/Paper (Model)/2 Revision 1/Figures New/Fig5r.png','-png','-transparent','-m2')
    close
end

%% Figure 6 new
if flags.do_fig6new
    figure('PaperUnits','centimeters','PaperType','A4','Paperposition',[0, 0, 21, 29.7],'Units','centimeters','Position',[0 0 21 29.7],'Resize','off')
    colormap(gray);
    xstepsize=2.5;
    ystepsize=5;
    methodLabel=['MAX';'CTD';'AGD';'MCM'];
    
    if flags.do_recalc
        data=data_ziegelwanger2014('SPHERE_ROT','recalc');
    else
        data=data_ziegelwanger2014('SPHERE_ROT');
    end
    for ii=1:length(data.results)
        p_onaxis{1}(:,:,ii)=data.results(ii).MAX{1}.p_onaxis;
        p_onaxis{2}(:,:,ii)=data.results(ii).CTD{1}.p_onaxis;
        p_onaxis{3}(:,:,ii)=data.results(ii).AGD{1}.p_onaxis;
        p_onaxis{4}(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
        p_onaxisbin{1}(:,1,ii)=data.results(ii).MAX{1}.p_onaxisbin;
        p_onaxisbin{2}(:,1,ii)=data.results(ii).CTD{1}.p_onaxisbin;
        p_onaxisbin{3}(:,1,ii)=data.results(ii).AGD{1}.p_onaxisbin;
        p_onaxisbin{4}(:,1,ii)=data.results(ii).MCM{1}.p_onaxisbin;
    end
    
    for ii=1:length(data.results)
        resnormS(ii,1)=mean([data.results(ii).MAX{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).MAX{1}.performance.on_axis{2}.resnormS]);
        resnormS(ii,2)=mean([data.results(ii).CTD{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).CTD{1}.performance.on_axis{2}.resnormS]);
        resnormS(ii,3)=mean([data.results(ii).AGD{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).AGD{1}.performance.on_axis{2}.resnormS]);
        resnormS(ii,4)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).MCM{1}.performance.on_axis{2}.resnormS]);
    end
    
    idx=find(abs(data.phi)~=40);
    
    fprintf(['Resnorm (MAX): ' num2str(mean(resnormS(:,1))) ' samples\n']);
    fprintf(['Resnorm (CTD): ' num2str(mean(resnormS(:,2))) ' samples\n']);
%     fprintf(['Resnorm (AGD): ' num2str(mean(resnormS(:,3))) ' samples\n']);
    fprintf(['Resnorm (MCM): ' num2str(mean(resnormS(:,4))) ' samples\n']);
    
    % radii
    varl=[[squeeze(p_onaxis{1}(1,1,idx)) ...
        squeeze(p_onaxis{2}(1,1,idx))...
        squeeze(p_onaxis{3}(1,1,idx)) ...
        squeeze(p_onaxis{4}(1,1,idx))]*1000 ...
        data.radius(idx)];
    varr=[[squeeze(p_onaxis{1}(1,2,idx)) ...
        squeeze(p_onaxis{2}(1,2,idx))...
        squeeze(p_onaxis{3}(1,2,idx)) ...
        squeeze(p_onaxis{4}(1,2,idx))]*1000 ...
        data.radius(idx)];
    var=[[squeeze(p_onaxisbin{1}(1,1,idx)) ...
        squeeze(p_onaxisbin{2}(1,1,idx))...
        squeeze(p_onaxisbin{3}(1,1,idx)) ...
        squeeze(p_onaxisbin{4}(1,1,idx))]*1000 ...
        data.radius(idx)];
    err=abs([varl(:,1:4); varr(:,1:4)]-[varl(:,5); varr(:,5)]*[1 1 1 1]);
    fprintf(['Radius (MAX): average err is ' num2str(mean(err(:,1))) ' cm\n']);
    fprintf(['              standard deviation is ' num2str(std(err(:,1))) ' cm\n']);
    fprintf(['Radius (CTD): average err is ' num2str(mean(err(:,2))) ' cm\n']);
    fprintf(['              standard deviation is ' num2str(std(err(:,2))) ' cm\n']);
%     fprintf(['Radius (AGD): average err is ' num2str(mean(err(:,3))) ' cm\n']);
%     fprintf(['              standard deviation is ' num2str(std(err(:,3))) ' cm\n']);
    fprintf(['Radius (MCM): average err is ' num2str(mean(err(:,4))) ' cm\n']);
    fprintf(['              standard deviation is ' num2str(std(err(:,4))) ' cm\n']);

    h(1)=subplot(631);
    hist(varl(:,[1 2 4]),0:xstepsize:200)
    xlim([70,105])
    xlabel('r in mm (left ear)')
    ylabel('Relative Frequency')
    set(h(1),'XTick',[77.5 87.5 97.5])
    h(2)=subplot(632);
    hist(varr(:,[1 2 4]),0:xstepsize:200)
    xlim([70,105])
    xlabel('r in mm (right ear)')
    set(h(2),'XTick',[77.5 87.5 97.5])
%     h(3)=subplot(633);
%     hist(var(:,[1 2 4]),0:xstepsize:200)
%     xlim([70,105])
%     xlabel('r in mm')
%     set(h(3),'XTick',[77.5 87.5 97.5])
    clear var varl varr

    %phi
    varl=[[squeeze(p_onaxis{1}(2,1,idx)) ...
        squeeze(p_onaxis{2}(2,1,idx))...
        squeeze(p_onaxis{3}(2,1,idx)) ...
        squeeze(p_onaxis{4}(2,1,idx))]/pi*180 ...
        -data.phi(idx)+ones(length(data.phi(idx)),1)*90];
    varr=[mod([squeeze(p_onaxis{1}(2,2,idx)) ...
        squeeze(p_onaxis{2}(2,2,idx))...
        squeeze(p_onaxis{3}(2,2,idx)) ...
        squeeze(p_onaxis{4}(2,2,idx))]/pi*180,360) ...
        mod(-data.phi(idx)-ones(length(data.phi(idx)),1)*90,360)];
    err=abs([varl(:,1:4); varr(:,1:4)]-[varl(:,5); varr(:,5)]*[1 1 1 1]);
    fprintf(['Phi (MAX): average err is ' num2str(mean(err(:,1))) 'deg\n']);
    fprintf(['           standard deviation is ' num2str(std(err(:,1))) 'deg\n']);
    fprintf(['Phi (CTD): average err is ' num2str(mean(err(:,2))) 'deg\n']);
    fprintf(['           standard deviation is ' num2str(std(err(:,2))) 'deg\n']);
%     fprintf(['Phi (AGD): average err is ' num2str(mean(err(:,3))) 'deg\n']);
%     fprintf(['           standard deviation is ' num2str(std(err(:,3))) 'deg\n']);
    fprintf(['Phi (MCM): average err is ' num2str(mean(err(:,4))) 'deg\n']);
    fprintf(['           standard deviation is ' num2str(std(err(:,4))) 'deg\n']);

    h(4)=subplot(6,3,4);
    hist(varl(:,[1 2 4]),50:ystepsize:130)
    xlim([55,125])
    xlabel('phi_e in degree (left ear)')
    ylabel('Relative Frequency')
    set(h(4),'XTick',60:10:120)
    h(5)=subplot(6,3,5);
    hist(varr(:,[1 2 4]),230:ystepsize:310)
    xlim([235,305])
    xlabel('phi_e in degree (right ear)')
    set(h(5),'XTick',240:10:300)
    clear varl varr

    %theta
    varl=[[squeeze(p_onaxis{1}(3,1,idx)) ...
        squeeze(p_onaxis{2}(3,1,idx))...
        squeeze(p_onaxis{3}(3,1,idx)) ...
        squeeze(p_onaxis{4}(3,1,idx))]/pi*180 ...
        data.theta(idx)];
    varr=[[squeeze(p_onaxis{1}(3,2,idx)) ...
        squeeze(p_onaxis{2}(3,2,idx))...
        squeeze(p_onaxis{3}(3,2,idx)) ...
        squeeze(p_onaxis{4}(3,2,idx))]/pi*180 ...
        -data.theta(idx)];
    err=abs([varl(:,1:4); varr(:,1:4)]-[varl(:,5); varr(:,5)]*[1 1 1 1]);
    fprintf(['Theta (MAX): average err is ' num2str(mean(err(:,1))) 'deg\n']);
    fprintf(['             standard deviation is ' num2str(std(err(:,1))) 'deg\n']);
    fprintf(['Theta (CTD): average err is ' num2str(mean(err(:,2))) 'deg\n']);
    fprintf(['             standard deviation is ' num2str(std(err(:,2))) 'deg\n']);
%     fprintf(['Theta (AGD): average err is ' num2str(mean(err(:,3))) 'deg\n']);
%     fprintf(['             standard deviation is ' num2str(std(err(:,3))) 'deg\n']);
    fprintf(['Theta (MCM): average err is ' num2str(mean(err(:,4))) 'deg\n']);
    fprintf(['             standard deviation is ' num2str(std(err(:,4))) 'deg\n']);
    
    h(6)=subplot(6,3,7);
    hist(varl(:,[1 2 4])-1,-20:xstepsize:20)
    xlim([-15,15])
    xlabel('\theta_e in degree (left ear)')
    ylabel('Relative Frequency')
    set(h(6),'XTick',-10:5:10)
    h(7)=subplot(6,3,8);
    hist(varr(:,[1 2 4])-1,-20:xstepsize:20)
    xlim([-15,15])
    xlabel('\theta_e in degree (right ear)')
    set(h(7),'XTick',-10:5:10)
    leg=legend('MAX','CTD','MCM');
    XYData=get(leg,'Position');
    XYData(1)=XYData(1)+XYData(3)/4*3;
    XYData(3)=XYData(3)/4;
    set(leg,'Position',XYData);
    leg1=findobj(leg,'type','text');
    set(leg1,'FontSize',8)
    clear varl varr leg
        
    for ii=[1:2 4:7]
        set(h(ii),'Linewidth',2)
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end
    
    export_fig('~/Dropbox/HRTF-TOA (1)/Publications/Paper (Model)/2 Revision 1/Figures New/Fig6.png','-png','-r300')
    close
end

%% Figure 7, Figure 12 new
if flags.do_fig7new || flags.do_fig12new
    
    hrtf={'ARI','CIPIC','LISTEN'};

    %-------------------------------Load Data----------------------------------
    for kk=1:length(hrtf)
        if flags.do_recalc
            data=data_ziegelwanger2014(hrtf{kk},'recalc');
        else
            data=data_ziegelwanger2014(hrtf{kk});
        end
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        temp1=zeros(size(data.results(1).MCM{1}.p_onaxis,1),size(data.results(1).MCM{1}.p_onaxis,2),length(data.results));
        temp3=zeros(size(data.results(1).MCM{1}.p_offaxis,1),size(data.results(1).MCM{1}.p_offaxis,2),length(data.results));
        for ii=1:length(data.results)
            temp1(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
            temp3(:,1:size(data.results(ii).MCM{1}.p_offaxis,2),ii)=data.results(ii).MCM{1}.p_offaxis;
            temp5(ii)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormS ...
                data.results(ii).MCM{1}.performance.on_axis{2}.resnormS]);
            temp7(ii)=mean([data.results(ii).MCM{1}.performance.off_axis{1}.resnormS ...
                data.results(ii).MCM{1}.performance.off_axis{2}.resnormS]);
            temp9(ii)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormP ...
                data.results(ii).MCM{1}.performance.on_axis{2}.resnormP]);
            temp11(ii)=mean([data.results(ii).MCM{1}.performance.off_axis{1}.resnormP ...
                data.results(ii).MCM{1}.performance.off_axis{2}.resnormP]);
            temp13(ii)=data.results(ii).MCM{1}.performance.outlierRate;
        end
        p_onaxis{kk,1}=temp1;
        p_offaxis{kk,1}=temp3;
        resnormS_onaxis{kk,1}=temp5;
        resnormS_offaxis{kk,1}=temp7;
        resnormP_onaxis{kk,1}=temp9;
        resnormP_offaxis{kk,1}=temp11;
        outlierRate{kk,2}=temp13;
    end
    
    for kk=1:length(hrtf)
        if flags.do_recalc
            data=data_ziegelwanger2014(hrtf{kk},'recalc');
        else
            data=data_ziegelwanger2014(hrtf{kk});
        end
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        temp1=zeros(size(data.results(1).MCM{2}.p_onaxis,1),size(data.results(1).MCM{2}.p_onaxis,2),length(data.results));
        temp3=zeros(size(data.results(1).MCM{2}.p_offaxis,1),size(data.results(1).MCM{2}.p_offaxis,2),length(data.results));
        for ii=1:length(data.results)
            temp1(:,:,ii)=data.results(ii).MCM{2}.p_onaxis;
            temp3(:,1:size(data.results(ii).MCM{2}.p_offaxis,2),ii)=data.results(ii).MCM{2}.p_offaxis;
            temp5(ii)=mean([data.results(ii).MCM{2}.performance.on_axis{1}.resnormS ...
                data.results(ii).MCM{2}.performance.on_axis{2}.resnormS]);
            temp7(ii)=mean([data.results(ii).MCM{2}.performance.off_axis{1}.resnormS ...
                data.results(ii).MCM{2}.performance.off_axis{2}.resnormS]);
            temp9(ii)=mean([data.results(ii).MCM{2}.performance.on_axis{1}.resnormP ...
                data.results(ii).MCM{2}.performance.on_axis{2}.resnormP]);
            temp11(ii)=mean([data.results(ii).MCM{2}.performance.off_axis{1}.resnormP ...
                data.results(ii).MCM{2}.performance.off_axis{2}.resnormP]);
            temp13(ii)=data.results(ii).MCM{2}.performance.outlierRate;
        end
        p_onaxis{kk,2}=temp1;
        p_offaxis{kk,2}=temp3;
        resnormS_onaxis{kk,2}=temp5;
        resnormS_offaxis{kk,2}=temp7;
        resnormP_onaxis{kk,2}=temp9;
        resnormP_offaxis{kk,2}=temp11;
        outlierRate{kk,2}=temp13;
    end
    
    if flags.do_fig7new %Figure 7
        
        figure('PaperUnits','centimeters','PaperType','A4','Paperposition',[0, 0, 21, 29.7],'Units','centimeters','Position',[0 0 21 29.7],'Resize','off')
        colormap(gray);
        xstepsize=2;
        ystepsize=2;
% % %         fprintf('\nOn-AxisBin Model:\n')
% % %         fprintf(['\nResnormS (all,without outlier-detection):      ' num2str(mean([resnormS_onaxisbin{1,1} resnormS_onaxisbin{2,1} resnormS_onaxisbin{3,1}])) ' samples\n']);
% % %         fprintf(['ResnormS (all,with outlier-detection:          ' num2str(mean([resnormS_onaxisbin{1,2} resnormS_onaxisbin{2,2} resnormS_onaxisbin{3,2}])) ' samples\n']);
% % %         fprintf(['ResnormP (all,with outlier-detection:          ' num2str(mean([resnormP_onaxisbin{1,2} resnormP_onaxisbin{2,2} resnormP_onaxisbin{3,2}])) ' samples\n']);
        
        fprintf('\nOn-Axis Model:\n')
        fprintf(['\nResnormS (all,without outlier-detection):      ' num2str(mean([resnormS_onaxis{1,1} resnormS_onaxis{2,1} resnormS_onaxis{3,1}])) ' samples\n']);
        fprintf(['ResnormS (all,with outlier-detection:          ' num2str(mean([resnormS_onaxis{1,2} resnormS_onaxis{2,2} resnormS_onaxis{3,2}])) ' samples\n']);
        fprintf(['ResnormP (all,with outlier-detection:          ' num2str(mean([resnormP_onaxis{1,2} resnormP_onaxis{2,2} resnormP_onaxis{3,2}])) ' samples\n']);
        
        fprintf(['\nResnormS (ARI,without outlier-detection):      ' num2str(mean(resnormS_onaxis{1,1})) ' samples\n']);
        fprintf(['ResnormS (ARI,with outlier-detection:          ' num2str(mean(resnormS_onaxis{1,2})) ' samples\n']);
        fprintf(['ResnormP (ARI,with outlier-detection:          ' num2str(mean(resnormP_onaxis{1,2})) ' samples\n']);
        
        fprintf(['\nResnormS (CIPIC,without outlier-detection):    ' num2str(mean(resnormS_onaxis{2,1})) ' samples\n']);
        fprintf(['ResnormS (CIPIC,with outlier-detection:        ' num2str(mean(resnormS_onaxis{2,2})) ' samples\n']);
        fprintf(['ResnormP (CIPIC,without outlier-detection):    ' num2str(mean(resnormP_onaxis{2,1})) ' samples\n']);
        fprintf(['ResnormP (CIPIC,with outlier-detection:        ' num2str(mean(resnormP_onaxis{2,2})) ' samples\n']);
        
        fprintf(['\nResnormS (LISTEN,without outlier-detection):   ' num2str(mean(resnormS_onaxis{3,1})) ' samples\n']);
        fprintf(['ResnormS (LISTEN,with outlier-detection:       ' num2str(mean(resnormS_onaxis{3,2})) ' samples\n']);
        fprintf(['ResnormP (LISTEN,with outlier-detection:       ' num2str(mean(resnormP_onaxis{3,2})) ' samples\n']);
        
        % radii
        fprintf(['\nr  (Avg): ' ...
            num2str(mean([mean(squeeze(p_onaxis{1,1}(1,:,:)*1000)) ...
            mean(squeeze(p_onaxis{2,1}(1,:,:)*1000)) ...
            mean(squeeze(p_onaxis{3,1}(1,:,:)*1000))])) ' mm\n'])
        fprintf(['r  (Std): ' ...
            num2str(std([reshape(p_onaxis{1,1}(1,:,:)*1000,1,numel(p_onaxis{1,1}(1,:,:))) ...
            reshape(p_onaxis{2,1}(1,:,:)*1000,1,numel(p_onaxis{2,1}(1,:,:))) ...
            reshape(p_onaxis{3,1}(1,:,:)*1000,1,numel(p_onaxis{3,1}(1,:,:)))])) ' mm\n'])
        fprintf(['r  Difference (Avg): ' ...
            num2str(mean([mean(abs(diff(squeeze(p_onaxis{1,1}(1,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_onaxis{2,1}(1,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_onaxis{3,1}(1,:,:)*1000))))])) ' mm\n'])
        fprintf(['r  Difference (Std): ' ...
            num2str(std([abs(diff(squeeze(p_onaxis{1,1}(1,:,:)*1000))) ...
            abs(diff(squeeze(p_onaxis{2,1}(1,:,:)*1000))) ...
            abs(diff(squeeze(p_onaxis{3,1}(1,:,:)*1000)))])) ' mm\n'])
        fprintf(['r  Difference (Max): ' ...
            num2str(max([max(abs(diff(squeeze(p_onaxis{1,1}(1,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_onaxis{2,1}(1,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_onaxis{3,1}(1,:,:)*1000))))])) ' mm\n'])
        temp=1;
        for kk=1:length(hrtf)
            var(temp:temp+size(p_onaxis{kk,1},3)-1,:)= ...
                squeeze(mean(p_onaxis{kk,1}(1,:,:)*1000,2));
            varl(temp:temp+size(p_onaxis{kk,1},3)-1,:)= ...
                squeeze(p_onaxis{kk,1}(1,1,:)*1000);
            varr(temp:temp+size(p_onaxis{kk,1},3)-1,:)= ...
                squeeze(p_onaxis{kk,1}(1,2,:)*1000);
            temp=size(var,1)+1;
        end
        
        h(1)=subplot(631);
        hist([varl(:,1) varr(:,1)],0:xstepsize:200)
        xlim([50,130])
        xlabel('r in mm')
        ylabel('Relative Frequency')
        h(2)=subplot(632);
        hist(abs(varr(:,1)-varl(:,1)),0:ystepsize:100)
        xlim([0,60])
        xlabel('\Delta r in mm')
        clear var varl varr;
        
        % phi_e
        fprintf(['\nphi_e left  (Avg): ' ...
            num2str(mean([mean(squeeze(p_onaxis{1,1}(2,1,:)*180/pi)) ...
            mean(squeeze(p_onaxis{2,1}(2,1,:)*180/pi)) ...
            mean(squeeze(p_onaxis{3,1}(2,1,:)*180/pi))])) '°\n'])
        fprintf(['phi_e right (Avg): ' ...
            num2str(mean([mean(squeeze(p_onaxis{1,1}(2,2,:)*180/pi)) ...
            mean(squeeze(p_onaxis{2,1}(2,2,:)*180/pi)) ...
            mean(squeeze(p_onaxis{3,1}(2,2,:)*180/pi))])) '°\n'])
        fprintf(['phi_e left  (Std): ' ...
            num2str(std([reshape(p_onaxis{1,1}(2,1,:)*180/pi,1,numel(p_onaxis{1,1}(2,1,:))) ...
            reshape(p_onaxis{2,1}(2,1,:)*180/pi,1,numel(p_onaxis{2,1}(2,1,:))) ...
            reshape(p_onaxis{3,1}(2,1,:)*180/pi,1,numel(p_onaxis{3,1}(2,1,:)))])) '°\n'])
        fprintf(['phi_e right (Std): ' ...
            num2str(std([reshape(p_onaxis{1,1}(2,2,:)*180/pi,1,numel(p_onaxis{1,1}(2,2,:))) ...
            reshape(p_onaxis{2,1}(2,2,:)*180/pi,1,numel(p_onaxis{2,1}(2,2,:))) ...
            reshape(p_onaxis{3,1}(2,2,:)*180/pi,1,numel(p_onaxis{3,1}(2,2,:)))])) '°\n'])
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,1},3)-1,:)= ...
                squeeze(p_onaxis{kk,1}(2,1,:))*180/pi;
            varr(temp:temp+size(p_onaxis{kk,1},3)-1,:)= ...
                mod(squeeze(p_onaxis{kk,1}(2,2,:))*180/pi,360);
            temp=size(varl,1)+1;
        end
        
        h(3)=subplot(634);
        hist(varl(:,1),45:xstepsize:135)
        xlim([45,135])
        xlabel('phi_e in degree (left ear)')
        ylabel('Relative Frequency')
        h(4)=subplot(635);
        hist(varr(:,1),225:xstepsize:315)
        xlim([225,315])
        xlabel('phi_e in degree (right ear)')
        clear varl varr
        
        % theta_e
        fprintf(['\ntheta_e left  (Avg): ' ...
            num2str(mean([mean(squeeze(p_onaxis{1,1}(3,1,:)*180/pi)) ...
            mean(squeeze(p_onaxis{2,1}(3,1,:)*180/pi)) ...
            mean(squeeze(p_onaxis{3,1}(3,1,:)*180/pi))])) '°\n'])
        fprintf(['theta_e right (Avg): ' ...
            num2str(mean([mean(squeeze(p_onaxis{1,1}(3,2,:)*180/pi)) ...
            mean(squeeze(p_onaxis{2,1}(3,2,:)*180/pi)) ...
            mean(squeeze(p_onaxis{3,1}(3,2,:)*180/pi))])) '°\n'])
        fprintf(['theta_e left  (Std): ' ...
            num2str(std([reshape(p_onaxis{1,1}(3,1,:)*180/pi,1,numel(p_onaxis{1,1}(3,1,:))) ...
            reshape(p_onaxis{2,1}(3,1,:)*180/pi,1,numel(p_onaxis{2,1}(3,1,:))) ...
            reshape(p_onaxis{3,1}(3,1,:)*180/pi,1,numel(p_onaxis{3,1}(3,1,:)))])) '°\n'])
        fprintf(['theta_e right (Std): ' ...
            num2str(std([reshape(p_onaxis{1,1}(3,2,:)*180/pi,1,numel(p_onaxis{1,1}(3,2,:))) ...
            reshape(p_onaxis{2,1}(3,2,:)*180/pi,1,numel(p_onaxis{2,1}(3,2,:))) ...
            reshape(p_onaxis{3,1}(3,2,:)*180/pi,1,numel(p_onaxis{3,1}(3,2,:)))])) '°\n\n'])
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,1},3)-1,:)= ...
                squeeze(p_onaxis{kk,1}(3,1,:))*180/pi;
            varr(temp:temp+size(p_onaxis{kk,1},3)-1,:)= ...
                squeeze(p_onaxis{kk,1}(3,2,:))*180/pi;
            temp=size(varl,1)+1;
        end
        
        h(5)=subplot(637);
        hist(varl(:,1),-45:xstepsize:45)
        xlim([-45,45])
        xlabel('\theta_e in degree (left ear)')
        ylabel('Relative Frequency')
        h(6)=subplot(638);
        hist(varr(:,1),-45:xstepsize:45)
        xlim([-45,45])
        xlabel('\theta_e in degree (right ear)')
        clear varl varr
        
        for ii=1:length(h)
            set(h(ii),'Linewidth',2)
            tmp=findobj(h(ii),'Type','patch');
            set(tmp,'EdgeColor','k');
        end
        
        export_fig('~/Dropbox/HRTF-TOA (1)/Publications/Paper (Model)/2 Revision 1/Figures New/Fig7.png','-png','-r300')
        close
    end
    
    if flags.do_fig12new %Figure12
        figure('PaperUnits','centimeters','PaperType','A4','Paperposition',[0, 0, 21, 29.7],'Units','centimeters','Position',[0 0 21 29.7],'Resize','off')
        colormap(gray);
        xstepsize=2;
        ystepsize=2;
        
% % %         fprintf('\nOff-AxisBin Model:\n')
% % %         fprintf(['\nResnormS (all,without outlier-detection):    ' num2str(mean([resnormS_offaxisbin{1,1} resnormS_offaxisbin{2,1} resnormS_offaxisbin{3,1}])) ' samples\n']);
% % %         fprintf(['ResnormS (all,with outlier-detection:        ' num2str(mean([resnormS_offaxisbin{1,2} resnormS_offaxisbin{2,2} resnormS_offaxisbin{3,2}])) ' samples\n']);
% % %         fprintf(['ResnormP (all,with outlier-detection:        ' num2str(mean([resnormP_offaxisbin{1,2} resnormP_offaxisbin{2,2} resnormP_offaxisbin{3,2}])) ' samples\n']);
        
        fprintf('\nOff-Axis Model:\n')
        fprintf(['\nResnormS (all,without outlier-detection):    ' num2str(mean([resnormS_offaxis{1,1} resnormS_offaxis{2,1} resnormS_offaxis{3,1}])) ' samples\n']);
        fprintf(['ResnormS (all,with outlier-detection:        ' num2str(mean([resnormS_offaxis{1,2} resnormS_offaxis{2,2} resnormS_offaxis{3,2}])) ' samples\n']);
        fprintf(['ResnormP (all,with outlier-detection:        ' num2str(mean([resnormP_offaxis{1,2} resnormP_offaxis{2,2} resnormP_offaxis{3,2}])) ' samples\n']);
        fprintf(['OutlierRate (all):                           ' num2str(mean([outlierRate{1,2} outlierRate{2,2} outlierRate{3,2}])) ' samples\n']);
        
        fprintf(['\nResnormS (ARI,without outlier-detection):    ' num2str(mean(resnormS_offaxis{1,1})) ' samples\n']);
        fprintf(['ResnormS (ARI,with outlier-detection:        ' num2str(mean(resnormS_offaxis{1,2})) ' samples\n']);
        fprintf(['ResnormP (ARI,with outlier-detection:        ' num2str(mean(resnormP_offaxis{1,2})) ' samples\n']);
        fprintf(['OutlierRate (ARI):                           ' num2str(mean(outlierRate{1,2})) ' samples\n']);
        
        fprintf(['\nResnormS (CIPIC,without outlier-detection):  ' num2str(mean(resnormS_offaxis{2,1})) ' samples\n']);
        fprintf(['ResnormS (CIPIC,with outlier-detection:      ' num2str(mean(resnormS_offaxis{2,2})) ' samples\n']);
        fprintf(['ResnormP (CIPIC,with outlier-detection:      ' num2str(mean(resnormP_offaxis{2,2})) ' samples\n']);
        fprintf(['OutlierRate (CIPIC):                         ' num2str(mean(outlierRate{2,2})) ' samples\n']);
        
        fprintf(['\nResnormS (LISTEN,without outlier-detection): ' num2str(mean(resnormS_offaxis{3,1})) ' samples\n']);
        fprintf(['ResnormS (LISTEN,with outlier-detection:     ' num2str(mean(resnormS_offaxis{3,2})) ' samples\n']);
        fprintf(['ResnormP (LISTEN,with outlier-detection:     ' num2str(mean(resnormP_offaxis{3,2})) ' samples\n']);
        fprintf(['OutlierRate (LISTEN):                        ' num2str(mean(outlierRate{3,2})) ' samples\n']);
        
        % radii
        fprintf(['\nr  (Avg):            ' ...
            num2str(mean([mean(squeeze(p_offaxis{1,2}(1,:,:)*1000)) ...
            mean(squeeze(p_offaxis{2,2}(1,:,:)*1000)) ...
            mean(squeeze(p_offaxis{3,2}(1,:,:)*1000))])) ' mm\n'])
        fprintf(['r  (Std):            ' ...
            num2str(std([reshape(p_offaxis{1,2}(1,:,:)*1000,1,numel(p_offaxis{1,2}(1,:,:))) ...
            reshape(p_offaxis{2,2}(1,:,:)*1000,1,numel(p_offaxis{2,2}(1,:,:))) ...
            reshape(p_offaxis{3,2}(1,:,:)*1000,1,numel(p_offaxis{3,2}(1,:,:)))])) ' mm\n'])
        fprintf(['r  Difference (Avg): ' ...
            num2str(mean([mean(abs(diff(squeeze(p_offaxis{1,2}(1,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{2,2}(1,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{3,2}(1,:,:)*1000))))])) ' mm\n'])
        fprintf(['r  Difference (Std): ' ...
            num2str(std([abs(diff(squeeze(p_offaxis{1,2}(1,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{2,2}(1,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{3,2}(1,:,:)*1000)))])) ' mm\n'])
        fprintf(['r  Difference (Max): ' ...
            num2str(max([max(abs(diff(squeeze(p_offaxis{1,2}(1,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{2,2}(1,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{3,2}(1,:,:)*1000))))])) ' mm\n'])
        
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(1,1,:)*1000);
            varr(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(1,2,:)*1000);
% % %             var(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
% % %                 squeeze(p_offaxisbin{kk,2}(1,:,:)*1000);
            temp=size(varl,1)+1;
        end
        
        h(1)=subplot(631);
        hist([varl(:,1) varr(:,1)],0:xstepsize:200)
        xlim([50,130])
        xlabel('r in mm')
        ylabel('Relative Frequency')
        h(2)=subplot(632);
        hist(abs(varr(:,1)-varl(:,1)),0:ystepsize:100)
        xlim([0,60])
        xlabel('\Delta r in mm')
        h(3)=subplot(633);
% % %         hist(var(:,1),0:xstepsize:200)
% % %         xlim([50,130])
% % %         xlabel('r in mm')
        clear var varl varr

        %lateral displacements
        % xM
        subplot(412)
        fprintf(['\nxM (Avg):            ' ...
            num2str(mean([mean(abs(squeeze(p_offaxis{1,2}(2,:,:)*1000))) ...
            mean(abs(squeeze(p_offaxis{2,2}(2,:,:)*1000))) ...
            mean(abs(squeeze(p_offaxis{3,2}(2,:,:)*1000)))])) ' mm\n'])
        fprintf(['xM (Std):            ' ...
            num2str(std([reshape(p_offaxis{1,2}(2,:,:)*1000,1,numel(p_offaxis{1,2}(2,:,:))) ...
            reshape(p_offaxis{2,2}(2,:,:)*1000,1,numel(p_offaxis{2,2}(2,:,:))) ...
            reshape(p_offaxis{3,2}(2,:,:)*1000,1,numel(p_offaxis{3,2}(2,:,:)))])) ' mm\n'])
        fprintf(['xM Difference (Avg): ' ...
            num2str(mean([mean(abs(diff(squeeze(p_offaxis{1,2}(2,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{2,2}(2,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{3,2}(2,:,:)*1000))))])) ' mm\n'])
        fprintf(['xM Difference (Std): ' ...
            num2str(std([abs(diff(squeeze(p_offaxis{1,2}(2,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{2,2}(2,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{3,2}(2,:,:)*1000)))])) ' mm\n'])
        fprintf(['xM Difference (Max): ' ...
            num2str(max([max(abs(diff(squeeze(p_offaxis{1,2}(2,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{2,2}(2,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{3,2}(2,:,:)*1000))))])) ' mm\n'])
             
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(2,1,:)*1000);
            varr(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(2,2,:)*1000);
% % %             var(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
% % %                 squeeze(p_offaxisbin{kk,2}(2,:,:)*1000);
            temp=size(varl,1)+1;
        end
        
        h(4)=subplot(634);
        hist([varl(:,1) varr(:,1)],-100:xstepsize:100)
        xlim([-40,40])
        xlabel('x_M in mm')
        ylabel('Relative Frequency')
        h(5)=subplot(635);
        hist(abs(varr(:,1)-varl(:,1)),0:ystepsize:100)
        xlim([0,60])
        xlabel('\Delta x_M in mm')
        h(6)=subplot(636);
% % %         hist(var,-100:xstepsize:100)
% % %         xlim([-40,40])
% % %         xlabel('x_M in mm')
        clear var varl varr
        
        % yM
        subplot(413)
        fprintf(['\nyM (Avg):            ' ...
            num2str(mean([mean(abs(squeeze(p_offaxis{1,2}(3,:,:)*1000))) ...
            mean(abs(squeeze(p_offaxis{2,2}(3,:,:)*1000))) ...
            mean(abs(squeeze(p_offaxis{3,2}(3,:,:)*1000)))])) ' mm\n'])
        fprintf(['yM (Std):            ' ...
            num2str(std([reshape(p_offaxis{1,2}(3,:,:)*1000,1,numel(p_offaxis{1,2}(3,:,:))) ...
            reshape(p_offaxis{2,2}(3,:,:)*1000,1,numel(p_offaxis{2,2}(3,:,:))) ...
            reshape(p_offaxis{3,2}(3,:,:)*1000,1,numel(p_offaxis{3,2}(3,:,:)))])) ' mm\n'])
        fprintf(['yM Difference (Avg): ' ...
            num2str(mean([mean(abs(diff(squeeze(p_offaxis{1,2}(3,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{2,2}(3,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{3,2}(3,:,:)*1000))))])) ' mm\n'])
        fprintf(['yM Difference (Std): ' ...
            num2str(std([abs(diff(squeeze(p_offaxis{1,2}(3,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{2,2}(3,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{3,2}(3,:,:)*1000)))])) ' mm\n'])
        fprintf(['yM Difference (Max): ' ...
            num2str(max([max(abs(diff(squeeze(p_offaxis{1,2}(3,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{2,2}(3,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{3,2}(3,:,:)*1000))))])) ' mm\n'])
        
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(3,1,:)*1000);
            varr(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(3,2,:)*1000);
% % %             var(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
% % %                 squeeze(p_offaxisbin{kk,2}(3,:,:)*1000);
            temp=size(varl,1)+1;
        end
        
        h(7)=subplot(637);
        hist([varl(:,1) varr(:,1)],-100:xstepsize:100)
        xlim([-40,40])
        xlabel('y_M in mm')
        ylabel('Relative Frequency')
        h(8)=subplot(638);
        hist(abs(varr(:,1)-varl(:,1)),0:ystepsize:100)
        xlim([0,60])
        xlabel('\Delta y_M in mm')
        h(9)=subplot(639);
% % %         hist(var,-100:xstepsize:100)
% % %         xlim([-40,40])
% % %         xlabel('y_M in mm')
        clear var varl varr

        % zM
        subplot(414)
        fprintf(['\nzM (Avg):            ' ...
            num2str(mean([mean(abs(squeeze(p_offaxis{1,2}(4,:,:)*1000))) ...
            mean(abs(squeeze(p_offaxis{2,2}(4,:,:)*1000))) ...
            mean(abs(squeeze(p_offaxis{3,2}(4,:,:)*1000)))])) ' mm\n'])
        fprintf(['zM (Std):            ' ...
            num2str(std([reshape(p_offaxis{1,2}(4,:,:)*1000,1,numel(p_offaxis{1,2}(4,:,:))) ...
            reshape(p_offaxis{2,2}(4,:,:)*1000,1,numel(p_offaxis{2,2}(4,:,:))) ...
            reshape(p_offaxis{3,2}(4,:,:)*1000,1,numel(p_offaxis{3,2}(4,:,:)))])) ' mm\n'])
        fprintf(['zM Difference (Avg): ' ...
            num2str(mean([mean(abs(diff(squeeze(p_offaxis{1,2}(4,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{2,2}(4,:,:)*1000)))) ...
            mean(abs(diff(squeeze(p_offaxis{3,2}(4,:,:)*1000))))])) ' mm\n'])
        fprintf(['zM Difference (Std): ' ...
            num2str(std([abs(diff(squeeze(p_offaxis{1,2}(4,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{2,2}(4,:,:)*1000))) ...
            abs(diff(squeeze(p_offaxis{3,2}(4,:,:)*1000)))])) ' mm\n'])
        fprintf(['zM Difference (Max): ' ...
            num2str(max([max(abs(diff(squeeze(p_offaxis{1,2}(4,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{2,2}(4,:,:)*1000)))) ...
            max(abs(diff(squeeze(p_offaxis{3,2}(4,:,:)*1000))))])) ' mm\n'])
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(4,1,:)*1000);
            varr(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(4,2,:)*1000);
% % %             var(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
% % %                 squeeze(p_offaxisbin{kk,2}(4,:,:)*1000);
            temp=size(varl,1)+1;
        end
        
        h(10)=subplot(6,3,10);
        hist([varl(:,1) varr(:,1)],-100:xstepsize:100)
        xlim([-40,40])
        xlabel('z_M in mm')
        ylabel('Relative Frequency')
        h(11)=subplot(6,3,11);
        hist(abs(varr(:,1)-varl(:,1)),0:ystepsize:100)
        xlim([0,60])
        xlabel('\Delta z_M in mm')
        h(12)=subplot(6,3,12);
% % %         hist(var,-100:xstepsize:100)
% % %         xlim([-40,40])
% % %         xlabel('z_M in mm')
        clear var varl varr

        %ear positions
        % phi_e
        fprintf(['\nphi_e left  (Avg):   ' ...
            num2str(mean([mean(squeeze(p_offaxis{1,1}(6,1,:)*180/pi)) ...
            mean(squeeze(p_offaxis{2,1}(6,1,:)*180/pi)) ...
            mean(squeeze(p_offaxis{3,1}(6,1,:)*180/pi))])) '°\n'])
        fprintf(['phi_e right (Avg):   ' ...
            num2str(mean([mean(squeeze(p_offaxis{1,1}(6,2,:)*180/pi)) ...
            mean(squeeze(p_offaxis{2,1}(6,2,:)*180/pi)) ...
            mean(squeeze(p_offaxis{3,1}(6,2,:)*180/pi))])) '°\n'])
        fprintf(['phi_e left  (Std):   ' ...
            num2str(std([reshape(p_offaxis{1,1}(6,1,:)*180/pi,1,numel(p_offaxis{1,1}(6,1,:))) ...
            reshape(p_offaxis{2,1}(6,1,:)*180/pi,1,numel(p_offaxis{2,1}(6,1,:))) ...
            reshape(p_offaxis{3,1}(6,1,:)*180/pi,1,numel(p_offaxis{3,1}(6,1,:)))])) '°\n'])
        fprintf(['phi_e right (Std):   ' ...
            num2str(std([reshape(p_offaxis{1,1}(6,2,:)*180/pi,1,numel(p_offaxis{1,1}(6,2,:))) ...
            reshape(p_offaxis{2,1}(6,2,:)*180/pi,1,numel(p_offaxis{2,1}(6,2,:))) ...
            reshape(p_offaxis{3,1}(6,2,:)*180/pi,1,numel(p_offaxis{3,1}(6,2,:)))])) '°\n'])
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(6,1,:))*180/pi;
            varr(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                mod(squeeze(p_offaxis{kk,2}(6,2,:))*180/pi,360);
            temp=size(varl,1)+1;
        end
        
        h(13)=subplot(6,3,13);
        hist(varl(:,1),45:xstepsize:135)
        xlim([45,135])
        xlabel('phi_e in degree (left ear)')
        ylabel('Relative Frequency')
        h(14)=subplot(6,3,14);
        hist(varr(:,1),225:xstepsize:315)
        xlim([225,315])
        xlabel('phi_e in degree (right ear)')
        clear varl varr

        % theta_e
        fprintf(['\nphi_e left  (Avg):   ' ...
            num2str(mean([mean(squeeze(p_offaxis{1,1}(7,1,:)*180/pi)) ...
            mean(squeeze(p_offaxis{2,1}(7,1,:)*180/pi)) ...
            mean(squeeze(p_offaxis{3,1}(7,1,:)*180/pi))])) '°\n'])
        fprintf(['phi_e right (Avg):   ' ...
            num2str(mean([mean(squeeze(p_offaxis{1,1}(7,2,:)*180/pi)) ...
            mean(squeeze(p_offaxis{2,1}(7,2,:)*180/pi)) ...
            mean(squeeze(p_offaxis{3,1}(7,2,:)*180/pi))])) '°\n'])
        fprintf(['phi_e left  (Std):   ' ...
            num2str(std([reshape(p_offaxis{1,1}(7,1,:)*180/pi,1,numel(p_offaxis{1,1}(7,1,:))) ...
            reshape(p_offaxis{2,1}(7,1,:)*180/pi,1,numel(p_offaxis{2,1}(6,1,:))) ...
            reshape(p_offaxis{3,1}(7,1,:)*180/pi,1,numel(p_offaxis{3,1}(6,1,:)))])) '°\n'])
        fprintf(['phi_e right (Std):   ' ...
            num2str(std([reshape(p_offaxis{1,1}(7,2,:)*180/pi,1,numel(p_offaxis{1,1}(7,2,:))) ...
            reshape(p_offaxis{2,1}(7,2,:)*180/pi,1,numel(p_offaxis{2,1}(6,2,:))) ...
            reshape(p_offaxis{3,1}(7,2,:)*180/pi,1,numel(p_offaxis{3,1}(6,2,:)))])) '°\n\n'])
        temp=1;
        for kk=1:length(hrtf)
            varl(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(7,1,:))*180/pi;
            varr(temp:temp+size(p_onaxis{kk,2},3)-1,:)= ...
                squeeze(p_offaxis{kk,2}(7,2,:))*180/pi;
            temp=size(varl,1)+1;
        end
        
        h(15)=subplot(6,3,16);
        hist(varl(:,1),-45:xstepsize:45)
        xlim([-45,45])
        xlabel('\theta_e in degree (left ear)')
        ylabel('Relative Frequency')
        h(16)=subplot(6,3,17);
        hist(varr(:,1),-45:xstepsize:45)
        xlim([-45,45])
        xlabel('\theta_e in degree (right ear)')
        clear varl varr
        
        for ii=1:length(h)
            set(h(ii),'Linewidth',2)
            tmp=findobj(h(ii),'Type','patch');
            set(tmp,'EdgeColor','k');
        end
        
        export_fig('~/Fig12.png','-png','-r300')
        export_fig('~/Dropbox/HRTF-TOA (1)/Publications/Paper (Model)/2 Revision 1/Figures New/Fig12.png','-png','-r300')
        close
    end
    
end

%% Figure 8 new
if flags.do_fig8new 
    figure('PaperUnits','centimeters','PaperType','A4','Paperposition',[0, 0, 21, 29.7],'Units','centimeters','Position',[0 0 21 29.7],'Resize','off')
    clear h
    colormap(gray);
    xstepsize=2.5;
    ystepsize=5;
    methodLabel=['MAX';'CTD';'AGD';'MCM'];
    
    if flags.do_recalc
        data=data_ziegelwanger2014('SPHERE_DIS','recalc');
    else
        data=data_ziegelwanger2014('SPHERE_DIS');
    end
    for ii=1:length(data.results)
        p_onaxis{1}(:,:,ii)=data.results(ii).MAX{1}.p_onaxis;
        p_onaxis{2}(:,:,ii)=data.results(ii).CTD{1}.p_onaxis;
        p_onaxis{3}(:,:,ii)=data.results(ii).AGD{1}.p_onaxis;
        p_onaxis{4}(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
    end
    
    % radii
    varl=[[squeeze(p_onaxis{1}(1,1,:)) ...
        squeeze(p_onaxis{2}(1,1,:))...
        squeeze(p_onaxis{3}(1,1,:)) ...
        squeeze(p_onaxis{4}(1,1,:))]*1000 ...
        data.radius(:)];
    varr=[[squeeze(p_onaxis{1}(1,2,:)) ...
        squeeze(p_onaxis{2}(1,2,:))...
        squeeze(p_onaxis{3}(1,2,:)) ...
        squeeze(p_onaxis{4}(1,2,:))]*1000 ...
        data.radius(:)];
    err=abs([varl(:,1:4); varr(:,1:4)]-[varl(:,5); varr(:,5)]*[1 1 1 1]);
    fprintf(['Radius: average err is ' num2str(mean(err(:,1))) ' cm\n']);
    fprintf(['        standard deviation is ' num2str(std(err(:,1))) ' cm\n']);

    h(1)=subplot(631);
    hist(varl(:,[1 2 4]),0:xstepsize:200)
    xlim([50,125])
    xlabel('r in mm (left ear)')
    set(gca,'XMinorTick','on')
    ylabel('Relative Frequency')
    leg=legend('MAX','CTD','MCM','Location','NorthWest');
    leg1=findobj(leg,'type','text');
    set(leg1,'FontSize',8)
    set(h(1),'XTick',[77.5 87.5 97.5])
    h(2)=subplot(632);
    hist(varr(:,[1 2 4]),0:xstepsize:200)
    xlim([50,125])
    xlabel('r in mm (right ear)')
    set(gca,'XMinorTick','on')
    set(h(2),'XTick',[77.5 87.5 97.5])
    clear var varl varr
    
    for ii=1:length(h)
        set(h(ii),'Linewidth',2)
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end
    
    export_fig('~/Dropbox/HRTF-TOA (1)/Publications/Paper (Model)/2 Revision 1/Figures New/Fig8.png','-png','-r300')
    close
end

%% Figure 10 new
if flags.do_fig10new
    figure('PaperUnits','centimeters','PaperType','A4','Paperposition',[0, 0, 21, 29.7],'Units','centimeters','Position',[0 0 21 29.7],'Resize','off')
    clear h
    colormap(gray);
    xstepsize=2.5;
    ystepsize=5;
    
    if flags.do_recalc
        data=data_ziegelwanger2014('SPHERE_DIS','recalc');
    else
        data=data_ziegelwanger2014('SPHERE_DIS');
    end
    for ii=1:length(data.results)
        p_offaxis{1}(:,:,ii)=data.results(ii).MAX{1}.p_offaxis;
        p_offaxis{2}(:,:,ii)=data.results(ii).CTD{1}.p_offaxis;
        p_offaxis{3}(:,:,ii)=data.results(ii).AGD{1}.p_offaxis;
        p_offaxis{4}(:,:,ii)=data.results(ii).MCM{1}.p_offaxis;
        p_offaxisbin{1}(:,:,ii)=data.results(ii).MAX{1}.p_offaxisbin(1:7);
        p_offaxisbin{2}(:,:,ii)=data.results(ii).CTD{1}.p_offaxisbin(1:7);
        p_offaxisbin{3}(:,:,ii)=data.results(ii).AGD{1}.p_offaxisbin(1:7);
        p_offaxisbin{4}(:,:,ii)=data.results(ii).MCM{1}.p_offaxisbin(1:7);
        resnormS_offaxis{1}(ii)=mean([data.results(ii).MAX{1}.performance.off_axis{1}.resnormS ...
            data.results(ii).MAX{1}.performance.off_axis{2}.resnormS]);
        resnormS_offaxisbin{1}(ii)=data.results(ii).MAX{1}.performance.off_axisbin.resnormS;
        resnormP_offaxis{1}(ii)=mean([data.results(ii).MAX{1}.performance.off_axis{1}.resnormP ...
            data.results(ii).MAX{1}.performance.off_axis{2}.resnormP]);
        resnormP_offaxisbin{1}(ii)=data.results(ii).MAX{1}.performance.off_axisbin.resnormP;
        resnormS_offaxis{2}(ii)=mean([data.results(ii).CTD{1}.performance.off_axis{1}.resnormS ...
            data.results(ii).CTD{1}.performance.off_axis{2}.resnormS]);
        resnormS_offaxisbin{2}(ii)=data.results(ii).CTD{1}.performance.off_axisbin.resnormS;
        resnormP_offaxis{2}(ii)=mean([data.results(ii).CTD{1}.performance.off_axis{1}.resnormP ...
            data.results(ii).CTD{1}.performance.off_axis{2}.resnormP]);
        resnormP_offaxisbin{2}(ii)=data.results(ii).CTD{1}.performance.off_axisbin.resnormP;
        resnormS_offaxis{3}(ii)=mean([data.results(ii).AGD{1}.performance.off_axis{1}.resnormS ...
            data.results(ii).AGD{1}.performance.off_axis{2}.resnormS]);
        resnormS_offaxisbin{3}(ii)=data.results(ii).AGD{1}.performance.off_axisbin.resnormS;
        resnormP_offaxis{3}(ii)=mean([data.results(ii).AGD{1}.performance.off_axis{1}.resnormP ...
            data.results(ii).AGD{1}.performance.off_axis{2}.resnormP]);
        resnormP_offaxisbin{3}(ii)=data.results(ii).AGD{1}.performance.off_axisbin.resnormP;
        resnormS_offaxis{4}(ii)=mean([data.results(ii).MCM{1}.performance.off_axis{1}.resnormS ...
            data.results(ii).MCM{1}.performance.off_axis{2}.resnormS]);
        resnormS_offaxisbin{4}(ii)=data.results(ii).MCM{1}.performance.off_axisbin.resnormS;
        resnormP_offaxis{4}(ii)=mean([data.results(ii).MCM{1}.performance.off_axis{1}.resnormP ...
            data.results(ii).MCM{1}.performance.off_axis{2}.resnormP]);
        resnormP_offaxisbin{4}(ii)=data.results(ii).MCM{1}.performance.off_axisbin.resnormP;
    end
    
    center=[data.xM(1:length(data.xM)/3) data.yM(1:length(data.yM)/3) data.zM(1:length(data.zM)/3)];
    [~,idx3]=sort(squeeze(center(:,3)));
    [~,idx2]=sort(squeeze(center(:,1)));
    [~,idx1]=sort(squeeze(center(:,2)));
    idx=idx3(idx2(idx1));
    idx=[idx; idx+length(data.xM)/3; idx+length(data.xM)/3*2];
    data.radius=data.radius(idx);
    data.xM=data.xM(idx);
    data.yM=data.yM(idx);
    data.zM=data.zM(idx);
    
    fprintf('\nOff-AxisBin Model:\n')
    fprintf(['\nResnormS (MAX): ' num2str(mean(resnormS_offaxisbin{1})) ' samples\n']);
    fprintf(['ResnormS (CTD): ' num2str(mean(resnormS_offaxisbin{2})) ' samples\n']);
    fprintf(['ResnormS (MCM): ' num2str(mean(resnormS_offaxisbin{4})) ' samples\n']);

    fprintf('\nOff-Axis Model:\n')
    fprintf(['\nResnormS (MAX): ' num2str(mean(resnormS_offaxis{1})) ' samples\n']);
    fprintf(['ResnormS (CTD): ' num2str(mean(resnormS_offaxis{2})) ' samples\n']);
    fprintf(['ResnormS (MCM): ' num2str(mean(resnormS_offaxis{4})) ' samples\n']);

    %radii
    varl=[[squeeze(p_offaxis{1}(1,1,:)) ...
        squeeze(p_offaxis{2}(1,1,:))...
        squeeze(p_offaxis{3}(1,1,:)) ...
        squeeze(p_offaxis{4}(1,1,:))]*1000 ...
        data.radius];
    varr=[[squeeze(p_offaxis{1}(1,2,:)) ...
        squeeze(p_offaxis{2}(1,2,:))...
        squeeze(p_offaxis{3}(1,2,:)) ...
        squeeze(p_offaxis{4}(1,2,:))]*1000 ...
        data.radius];
    var=[[squeeze(p_offaxisbin{1}(1,:,:)) ...
        squeeze(p_offaxisbin{2}(1,:,:))...
        squeeze(p_offaxisbin{3}(1,:,:)) ...
        squeeze(p_offaxisbin{4}(1,:,:))]*1000 ...
        data.radius];
%     err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
%     err=reshape(err,numel(err),1);
%     temp1=err;
%     fprintf(['Radius: average err is ' num2str(mean(err)) ' cm\n']);
%     fprintf(['        standard deviation is ' num2str(std(err)) ' cm\n']);
%     fprintf(['        maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
%     fprintf(['        average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);

    h(1)=subplot(631);
    hist([varl(:,[1 2 4]); varr(:,[1 2 4])],0:xstepsize:200)
    xlim([70,105])
    xlabel('r in mm')
    ylabel('Relative Frequency')
    set(h(1),'XTick',[77.5 87.5 97.5])
    h(2)=subplot(632);
    hist(abs(varr(:,[1 2 4])-varl(:,[1 2 4])),0:xstepsize:200)
    xlim([-5,35])
    xlabel('\Delta r in mm')
    legend('MAX','CTD','MCM')
    set(h(2),'XTick',[0 10 20 30])
    h(3)=subplot(633);
    hist(var(:,[1 2 4]),0:xstepsize:200)
    xlim([70,105])
    xlabel('r in mm')
    set(h(3),'XTick',[77.5 87.5 97.5])
    clear var varl varr
    

    %xM
    varl=[[squeeze(p_offaxis{1}(2,1,:)) ...
        squeeze(p_offaxis{2}(2,1,:))...
        squeeze(p_offaxis{3}(2,1,:)) ...
        squeeze(p_offaxis{4}(2,1,:))]*1000 ...
        -data.xM*1000];
    varr=[[squeeze(p_offaxis{1}(2,2,:)) ...
        squeeze(p_offaxis{2}(2,2,:))...
        squeeze(p_offaxis{3}(2,2,:)) ...
        squeeze(p_offaxis{4}(2,2,:))]*1000 ...
        -data.xM*1000];
    var=[[squeeze(p_offaxisbin{1}(2,:,:)) ...
        squeeze(p_offaxisbin{2}(2,:,:))...
        squeeze(p_offaxisbin{3}(2,:,:)) ...
        squeeze(p_offaxisbin{4}(2,:,:))]*1000 ...
        -data.xM*1000];
%     err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
%     err=reshape(err,numel(err),1);
%     fprintf(['xM: average err is ' num2str(mean(err)) ' cm\n']);
%     fprintf(['    standard deviation is ' num2str(std(err)) ' cm\n']);
%     fprintf(['    maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
%     fprintf(['    average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        
    h(4)=subplot(634);
    hist([varl(:,[1 2 4]); varr(:,[1 2 4])],-100:xstepsize:100)
    xlim([-30,10])
    xlabel('x_M in mm')
    ylabel('Relative Frequency')
    set(h(4),'XTick',[-20 -10 0])
    h(5)=subplot(635);
    hist(abs(varr(:,[1 2 4])-varl(:,[1 2 4])),0:xstepsize:100)
    xlim([-5,35])
    xlabel('\Delta x_M in mm')
    set(h(5),'XTick',[0 10 20 30])    
    h(6)=subplot(636);
    hist(var(:,[1 2 4]),-100:xstepsize:100)
    xlim([-30,10])
    xlabel('x_M in mm')
    set(h(6),'XTick',[-20 -10 0])
    clear var varl varr
    

    %yM
    varl=[[squeeze(p_offaxis{1}(3,1,:)) ...
        squeeze(p_offaxis{2}(3,1,:))...
        squeeze(p_offaxis{3}(3,1,:)) ...
        squeeze(p_offaxis{4}(3,1,:))]*1000 ...
        -data.yM*1000];
    varr=[[squeeze(p_offaxis{1}(3,2,:)) ...
        squeeze(p_offaxis{2}(3,2,:))...
        squeeze(p_offaxis{3}(3,2,:)) ...
        squeeze(p_offaxis{4}(3,2,:))]*1000 ...
        -data.yM*1000];
    var=[[squeeze(p_offaxisbin{1}(3,:,:)) ...
        squeeze(p_offaxisbin{2}(3,:,:))...
        squeeze(p_offaxisbin{3}(3,:,:)) ...
        squeeze(p_offaxisbin{4}(3,:,:))]*1000 ...
        -data.yM*1000];
%     err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
%     err=reshape(err,numel(err),1);
%     temp1=[temp1; err];
%     fprintf(['yM: average err is ' num2str(mean(err)) ' cm\n']);
%     fprintf(['    standard deviation is ' num2str(std(err)) ' cm\n']);
%     fprintf(['    maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
%     fprintf(['    average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        
    h(7)=subplot(637);
    hist([varl(:,[1 2 4]); varr(:,[1 2 4])],-100:xstepsize:100)
    xlim([-10,30])
    xlabel('y_M in mm')
    ylabel('Relative Frequency')
    set(h(7),'XTick',[0 10 20])
    h(8)=subplot(638);
    hist(abs(varr(:,[1 2 4])-varl(:,[1 2 4])),0:xstepsize:100)
    xlim([-5,35])
    xlabel('\Delta y_M in mm')
    set(h(8),'XTick',[0 10 20 30])    
    h(9)=subplot(639);
    hist(var(:,[1 2 4]),-100:xstepsize:100)
    xlim([-10,30])
    xlabel('y_M in mm')
    set(h(9),'XTick',[0 10 20])
    clear var varl varr
    

    %zM
    varl=[[squeeze(p_offaxis{1}(4,1,:)) ...
        squeeze(p_offaxis{2}(4,1,:))...
        squeeze(p_offaxis{3}(4,1,:)) ...
        squeeze(p_offaxis{4}(4,1,:))]*1000 ...
        -data.zM*1000];
    varr=[[squeeze(p_offaxis{1}(4,2,:)) ...
        squeeze(p_offaxis{2}(4,2,:))...
        squeeze(p_offaxis{3}(4,2,:)) ...
        squeeze(p_offaxis{4}(4,2,:))]*1000 ...
        -data.zM*1000];
    var=[[squeeze(p_offaxisbin{1}(4,:,:)) ...
        squeeze(p_offaxisbin{2}(4,:,:))...
        squeeze(p_offaxisbin{3}(4,:,:)) ...
        squeeze(p_offaxisbin{4}(4,:,:))]*1000 ...
        -data.zM*1000];
%     err=abs([var(:,1) var(:,2)]-[var(:,3) var(:,3)]);
%     err=reshape(err,numel(err),1);
%     temp1=[temp1; err];
%     fprintf(['zM: average err is ' num2str(mean(err)) ' cm\n']);
%     fprintf(['    standard deviation is ' num2str(std(err)) ' cm\n']);
%     fprintf(['    maximum is ' num2str(max(abs(var(:,1)-var(:,2)))) ' cm\n']);
%     fprintf(['    average is ' num2str(mean(abs(var(:,1)-var(:,2)))) ' cm\n']);
        
    h(10)=subplot(6,3,10);
    hist([varl(:,[1 2 4]); varr(:,[1 2 4])],-100:xstepsize:100)
    xlim([-20,10])
    xlabel('z_M in mm')
    ylabel('Relative Frequency')
    set(h(10),'XTick',[-10 0])
    h(11)=subplot(6,3,11);
    hist(abs(varr(:,[1 2 4])-varl(:,[1 2 4])),0:xstepsize:100)
    xlim([-5,35])
    xlabel('\Delta z_M in mm')
    set(h(11),'XTick',[0 10 20 30])    
    h(12)=subplot(6,3,12);
    hist(var(:,[1 2 4]),-100:xstepsize:100)
    xlim([-20,10])
    xlabel('z_M in mm')
    set(h(12),'XTick',[-10 0])
    clear var varl varr
    


%     fprintf(['offset: average err is ' num2str(mean(temp1)) ' cm\n']);
%     fprintf(['        standard deviation is ' num2str(std(temp1)) ' cm\n']);
%     fprintf(['        maximum is ' num2str(max(abs(temp1))) ' cm\n']);

    % phi_e
    varl=[squeeze(p_offaxis{1}(6,1,:)) ...
        squeeze(p_offaxis{2}(6,1,:))...
        squeeze(p_offaxis{3}(6,1,:)) ...
        squeeze(p_offaxis{4}(6,1,:))]*180/pi ...
        ;
    varr=mod([squeeze(p_offaxis{1}(6,2,:)) ...
        squeeze(p_offaxis{2}(6,2,:))...
        squeeze(p_offaxis{3}(6,2,:)) ...
        squeeze(p_offaxis{4}(6,2,:))]*180/pi,360) ...
        ;

    h(13)=subplot(6,3,13);
    hist(varl(:,[1 2 4]),45:ystepsize:135)
    xlim([60,120])
    xlabel('phi_e in degree (left ear)')
    ylabel('Relative Frequency')
    h(14)=subplot(6,3,14);
    hist(varr(:,[1 2 4]),225:ystepsize:315)
    xlim([240,300])
    xlabel('phi_e in degree (right ear)')
%         ylabel('Relative Frequency')
    clear varl varr

    % theta_e
    varl=[squeeze(p_offaxis{1}(7,1,:)) ...
        squeeze(p_offaxis{2}(7,1,:))...
        squeeze(p_offaxis{3}(7,1,:)) ...
        squeeze(p_offaxis{4}(7,1,:))]*180/pi ...
        ;
    varr=[squeeze(p_offaxis{1}(6,2,:)) ...
        squeeze(p_offaxis{2}(6,2,:))...
        squeeze(p_offaxis{3}(6,2,:)) ...
        squeeze(p_offaxis{4}(6,2,:))] ...
        ;

    h(15)=subplot(6,3,16);
    hist(varl(:,[1 2 4]),-45:ystepsize:45)
    xlim([-30,30])
    xlabel('\theta_e in degree (left ear)')
    ylabel('Relative Frequency')
    h(16)=subplot(6,3,17);
    hist(varr(:,[1 2 4]),-45:ystepsize:45)
    xlim([-30,30])
    xlabel('\theta_e in degree (right ear)')
%         ylabel('Relative Frequency')
    clear varl varr

    for ii=1:length(h)
        set(h(ii),'Linewidth',2)
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end

    export_fig('~/Dropbox/HRTF-TOA (1)/Publications/Paper (Model)/2 Revision 1/Figures New/Fig10.png','-png','-r300')
    close
end


%% Table 2 new
if flags.do_tab2
    cd '/Volumes/ARI/Simulations/HRTF-TOA/SphereTorsoPinna/'
    method=1;
    methodLabel=['MAX';'CTD';'AGD';'MCM'];

    %%-Sphere------------------------------------------------------------------
    load(['Sphere' filesep 'hrtf_M_hrtf_2ears.mat']);
    Obj1=SOFAconvertARI2SOFA(hM,meta,stimPar);
    pos=zeros(Obj1.DimSize.M,8);
    pos(:,1:2)=Obj1.SourcePosition(:,1:2);
    [pos(:,6),pos(:,7)]=sph2hor(Obj1.SourcePosition(:,1),Obj1.SourcePosition(:,2));
    clear hM meta stimPar
    fprintf(['\n\nSphere:\n']);
    if flags.do_tab2
        for method=1:4
            [Obj1,results]=ziegelwanger2014(Obj1,method,0,0);
            %--------------------------------------------------------------------------
            [~,idxHor]=sort(pos(:,1));
            epsilon=5;
            hor_slope=zeros(Obj1.DimSize.M,1);
            for ele=min(pos(:,2)):epsilon:max(pos(:,2)) %calculate slope for each elevation along azimuth
                idx=find(pos(idxHor,2)>ele-epsilon/2 & pos(idxHor,2)<=ele+epsilon/2);
                if numel(idx)>1
                    idx(length(idx)+1)=idx(1);
                    hor_slope(idxHor(idx(1:end-1)),1)=diff(Obj1.Data.Delay(idxHor(idx),1))./(abs(abs(abs(diff(pos(idxHor(idx),1)))-180)-180)+0.00000000001);
                end
            end
            hor_sloperms=sqrt(sum(hor_slope.^2)/length(hor_slope));
            %--------------------------------------------------------------------------
            epsilon=2;
            sag_dev=zeros(Obj1.DimSize.M,1);
            sag_mean=sag_dev;
            for lat=-90:epsilon:90
                idx=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2);
                if length(idx)>2
                    sag_mean(idx,1)=mean(Obj1.Data.Delay(idx,1));
                    sag_dev(idx,1)=Obj1.Data.Delay(idx,1)-mean(Obj1.Data.Delay(idx,1));
                end
            end
            sag_min=min(sag_dev);
            sag_max=max(sag_dev);
            sag_var=sqrt(sum(sag_dev.^2)/length(sag_dev));
            %--------------------------------------------------------------------------
            fprintf(['\n' methodLabel(method,:) '\n']);
            fprintf(['Slope in horizontal planes (RMS): ' num2str(hor_sloperms) ' samples\n']);
            fprintf(['Variance in sagittal planes: ' num2str(sag_var) ' samples\n']);
            fprintf(['Min deviation in sagittal planes: ' num2str(sag_min) ' samples\n']);
            fprintf(['Max deviation in sagittal planes: ' num2str(sag_max) ' samples\n']);
            %--------------------------------------------------------------------------
            [Obj1,results]=ziegelwanger2014(Obj1,method,0,1);
            fprintf(['Radius: ' num2str(results.p_onaxis(1,1)*1000) ', ' num2str(results.p_onaxis(1,2)*1000) ' mm\n'])
            fprintf(['phi_e: ' num2str(results.p_onaxis(2,1)*180/pi) ', ' num2str(results.p_onaxis(2,2)*180/pi) ' °\n'])
            fprintf(['theta_e: ' num2str(results.p_onaxis(3,1)*180/pi) ', ' num2str(results.p_onaxis(3,2)*180/pi) ' °\n'])
            fprintf(['Resnorm: ' num2str(results.performance.on_axis{1}.resnormS) ' samples\n']);
            fprintf(['Resnorm: ' num2str(results.performance.off_axis{1}.resnormS) ' samples\n']);
        end
    end



    %%-SphereTorso-------------------------------------------------------------
    load(['SphereTorso' filesep 'hrtf_M_hrtf_2ears.mat']);
    Obj2=SOFAconvertARI2SOFA(hM,meta,stimPar);
    clear hM meta stimPar
    fprintf(['\n\nSphereTorso:\n\n']);
    if flags.do_tab2
        for method=1:4
            [Obj2,results]=ziegelwanger2014(Obj2,method,0,0);
            %--------------------------------------------------------------------------
            [~,idxHor]=sort(pos(:,1));
            epsilon=5;
            hor_slope=zeros(Obj2.DimSize.M,1);
            for ele=min(pos(:,2)):epsilon:max(pos(:,2)) %calculate slope for each elevation along azimuth
                idx=find(pos(idxHor,2)>ele-epsilon/2 & pos(idxHor,2)<=ele+epsilon/2);
                if numel(idx)>1
                    idx(length(idx)+1)=idx(1);
                    hor_slope(idxHor(idx(1:end-1)),1)=diff(Obj2.Data.Delay(idxHor(idx),1))./(abs(abs(abs(diff(pos(idxHor(idx),1)))-180)-180)+0.00000000001);
                end
            end
            hor_sloperms=sqrt(sum(hor_slope.^2)/length(hor_slope));
            %--------------------------------------------------------------------------
            epsilon=2;
            sag_dev=zeros(Obj2.DimSize.M,1);
            sag_mean=sag_dev;
            for lat=-90:epsilon:90
                idx=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2);
                if length(idx)>2
                    sag_mean(idx,1)=mean(Obj2.Data.Delay(idx,1));
                    sag_dev(idx,1)=Obj2.Data.Delay(idx,1)-mean(Obj2.Data.Delay(idx,1));
                end
            end
            sag_min=min(sag_dev);
            sag_max=max(sag_dev);
            sag_var=sqrt(sum(sag_dev.^2)/length(sag_dev));
            %--------------------------------------------------------------------------
            fprintf(['\n' methodLabel(method,:) '\n']);
            fprintf(['Slope in horizontal planes (RMS): ' num2str(hor_sloperms) ' samples\n']);
            fprintf(['Variance in sagittal planes: ' num2str(sag_var) ' samples\n']);
            fprintf(['Min deviation in sagittal planes: ' num2str(sag_min) ' samples\n']);
            fprintf(['Max deviation in sagittal planes: ' num2str(sag_max) ' samples\n']);
            %--------------------------------------------------------------------------
            [Obj2,results]=ziegelwanger2014(Obj2,method,0,1);
            fprintf(['Radius: ' num2str(results.p_onaxis(1,1)*1000) ', ' num2str(results.p_onaxis(1,2)*1000) ' mm\n'])
            fprintf(['phi_e: ' num2str(results.p_onaxis(2,1)*180/pi) ', ' num2str(results.p_onaxis(2,2)*180/pi) ' °\n'])
            fprintf(['theta_e: ' num2str(results.p_onaxis(3,1)*180/pi) ', ' num2str(results.p_onaxis(3,2)*180/pi) ' °\n'])
            fprintf(['Resnorm: ' num2str(results.performance.on_axis{1}.resnormS) ' samples\n']);
            fprintf(['Resnorm: ' num2str(results.performance.off_axis{1}.resnormS) ' samples\n']);
        end
    end



    %%-SphereTorsoPinna--------------------------------------------------------
    load(['SphereTorsoPinna' filesep 'hrtf_M_hrtf_2ears.mat']);
    Obj3=SOFAconvertARI2SOFA(hM,meta,stimPar);
    clear hM meta stimPar
    fprintf(['\n\nSphereTorsoPinna:\n']);
    if flags.do_tab2
        for method=1:4
            [Obj3,results]=ziegelwanger2014(Obj3,method,0,0);
            %--------------------------------------------------------------------------
            [~,idxHor]=sort(pos(:,1));
            epsilon=5;
            hor_slope=zeros(Obj3.DimSize.M,1);
            for ele=min(pos(:,2)):epsilon:max(pos(:,2)) %calculate slope for each elevation along azimuth
                idx=find(pos(idxHor,2)>ele-epsilon/2 & pos(idxHor,2)<=ele+epsilon/2);
                if numel(idx)>1
                    idx(length(idx)+1)=idx(1);
                    hor_slope(idxHor(idx(1:end-1)),1)=diff(Obj3.Data.Delay(idxHor(idx),1))./(abs(abs(abs(diff(pos(idxHor(idx),1)))-180)-180)+0.00000000001);
                end
            end
            hor_sloperms=sqrt(sum(hor_slope.^2)/length(hor_slope));
            %--------------------------------------------------------------------------
            epsilon=2;
            sag_dev=zeros(Obj3.DimSize.M,1);
            sag_mean=sag_dev;
            for lat=-90:epsilon:90
                idx=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2);
                if length(idx)>2
                    sag_mean(idx,1)=mean(Obj3.Data.Delay(idx,1));
                    sag_dev(idx,1)=Obj3.Data.Delay(idx,1)-mean(Obj3.Data.Delay(idx,1));
                end
            end
            sag_min=min(sag_dev);
            sag_max=max(sag_dev);
            sag_var=sqrt(sum(sag_dev.^2)/length(sag_dev));
            %--------------------------------------------------------------------------
            fprintf(['\n' methodLabel(method,:) '\n']);
            fprintf(['Slope in horizontal planes (RMS): ' num2str(hor_sloperms) ' samples\n']);
            fprintf(['Variance in sagittal planes: ' num2str(sag_var) ' samples\n']);
            fprintf(['Min deviation in sagittal planes: ' num2str(sag_min) ' samples\n']);
            fprintf(['Max deviation in sagittal planes: ' num2str(sag_max) ' samples\n']);
            %--------------------------------------------------------------------------
            [Obj3,results]=ziegelwanger2014(Obj3,method,0,1);
            fprintf(['Radius: ' num2str(results.p_onaxis(1,1)*1000) ', ' num2str(results.p_onaxis(1,2)*1000) ' mm\n'])
            fprintf(['phi_e: ' num2str(results.p_onaxis(2,1)*180/pi) ', ' num2str(results.p_onaxis(2,2)*180/pi) ' °\n'])
            fprintf(['theta_e: ' num2str(results.p_onaxis(3,1)*180/pi) ', ' num2str(results.p_onaxis(3,2)*180/pi) ' °\n'])
            fprintf(['Resnorm: ' num2str(results.performance.on_axis{1}.resnormS) ' samples\n']);
            fprintf(['Resnorm: ' num2str(results.performance.off_axis{1}.resnormS) ' samples\n']);
        end
    end



    % % % %%-Simulation-------------------------------------------------------------
    % % % load(['Simulation' filesep 'hrtf_M_hrtf.mat']);
    % % % meta.lat=0;
    % % % Obj5=SOFAconvertARI2SOFA(hM,meta,stimPar);
    % % % clear hM meta stimPar
    % % % for method=1:4
    % % %     [Obj5,results]=ziegelwanger2014(Obj5,method,0,0);
    % % %     %--------------------------------------------------------------------------
    % % %     [~,idxHor]=sort(pos(:,1));
    % % %     epsilon=5;
    % % %     hor_slope=zeros(Obj5.DimSize.M,1);
    % % %     for ele=min(pos(:,2)):epsilon:max(pos(:,2)) %calculate slope for each elevation along azimuth
    % % %         idx=find(pos(idxHor,2)>ele-epsilon/2 & pos(idxHor,2)<=ele+epsilon/2);
    % % %         if numel(idx)>1
    % % %             idx(length(idx)+1)=idx(1);
    % % %             hor_slope(idxHor(idx(1:end-1)),1)=diff(Obj5.Data.Delay(idxHor(idx),1))./(abs(abs(abs(diff(pos(idxHor(idx),1)))-180)-180)+0.00000000001);
    % % %         end
    % % %     end
    % % %     hor_sloperms=sqrt(sum(hor_slope.^2)/length(hor_slope));
    % % %     %--------------------------------------------------------------------------
    % % %     epsilon=2;
    % % %     sag_dev=zeros(Obj5.DimSize.M,1);
    % % %     sag_mean=sag_dev;
    % % %     for lat=-90:epsilon:90
    % % %         idx=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2);
    % % %         if length(idx)>2
    % % %             sag_mean(idx,1)=mean(Obj5.Data.Delay(idx,1));
    % % %             sag_dev(idx,1)=Obj5.Data.Delay(idx,1)-mean(Obj5.Data.Delay(idx,1));
    % % %         end
    % % %     end
    % % %     sag_min=min(sag_dev);
    % % %     sag_max=max(sag_dev);
    % % %     sag_var=sqrt(sum(sag_dev.^2)/length(sag_dev));
    % % %     %--------------------------------------------------------------------------
    % % %     fprintf(['\n\nSimulation:\n\n']);
    % % %     fprintf(['Slope in horizontal planes (RMS): ' num2str(hor_sloperms) ' samples\n']);
    % % %     fprintf(['Variance in sagittal planes: ' num2str(sag_var) ' samples\n']);
    % % %     fprintf(['Min deviation in sagittal planes: ' num2str(sag_min) ' samples\n']);
    % % %     fprintf(['Max deviation in sagittal planes: ' num2str(sag_max) ' samples\n']);
    % % %     %--------------------------------------------------------------------------
    % % % end

    %%-Measurement-------------------------------------------------------------
    hrtf={'ARI','CIPIC','LISTEN'};

    %-------------------------------Load Data----------------------------------
    for kk=1:length(hrtf)
        if flags.do_recalc
            data=data_ziegelwanger2014(hrtf{kk},'recalc');
        else
            data=data_ziegelwanger2014(hrtf{kk});
        end
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        if kk==1
            results=data.results;
        else
            results=[results data.results];
        end
    end
    clear data

    %MAX-------------------------------------------------------------------
    for ii=1:length(results)
        sloperms(ii)=mean(results(ii).MAX{2}.performance.sloperms(:));
        sag_var(ii)=mean(results(ii).MAX{2}.performance.sag_var);
        sag_max(ii)=max(results(ii).MAX{2}.performance.sag_max);
        resnormPon(ii)=mean([results(ii).MAX{1}.performance.on_axis{1}.resnormP ...
            results(ii).MAX{1}.performance.on_axis{2}.resnormP]);
        resnormPoff(ii)=mean([results(ii).MAX{2}.performance.off_axis{1}.resnormP ...
            results(ii).MAX{2}.performance.off_axis{2}.resnormP]);
        resnormSon(ii)=mean([results(ii).MAX{1}.performance.on_axis{1}.resnormS ...
            results(ii).MAX{1}.performance.on_axis{2}.resnormS]);
        resnormSoff(ii)=mean([results(ii).MAX{2}.performance.off_axis{1}.resnormS ...
            results(ii).MAX{2}.performance.off_axis{2}.resnormS]);
    end
    %----------------------------------------------------------------------
    fprintf(['\n\nMeasurement:\n\n']);
    fprintf(['MAX:\n']);
    fprintf(['Slope in horizontal planes (RMS): ' num2str(mean(sloperms)) ' samples\n']);
    fprintf(['Variance in sagittal planes: ' num2str(mean(sag_var)) ' samples\n']);
    fprintf(['Max deviation in sagittal planes: ' num2str(max(sag_max)) ' samples\n']);
    fprintf(['Resnorm P (on-axis): ' num2str(mean(resnormPon)) ' samples\n']);
    fprintf(['Resnorm P (off-axis):' num2str(mean(resnormPoff)) ' samples\n']);
    fprintf(['Resnorm S (on-axis): ' num2str(mean(resnormSon)) ' samples\n']);
    fprintf(['Resnorm S (off-axis):' num2str(mean(resnormSoff)) ' samples\n']);

    %CTD-------------------------------------------------------------------
    for ii=1:length(results)
        sloperms(ii)=mean(results(ii).CTD{2}.performance.sloperms);
        sag_var(ii)=mean(results(ii).CTD{2}.performance.sag_var);
        sag_max(ii)=max(results(ii).CTD{2}.performance.sag_max);
        resnormPon(ii)=mean([results(ii).CTD{1}.performance.on_axis{1}.resnormP ...
            results(ii).CTD{1}.performance.on_axis{2}.resnormP]);
        resnormPoff(ii)=mean([results(ii).CTD{2}.performance.off_axis{1}.resnormP ...
            results(ii).CTD{2}.performance.off_axis{2}.resnormP]);
        resnormSon(ii)=mean([results(ii).CTD{1}.performance.on_axis{1}.resnormS ...
            results(ii).CTD{1}.performance.on_axis{2}.resnormS]);
        resnormSoff(ii)=mean([results(ii).CTD{2}.performance.off_axis{1}.resnormS ...
            results(ii).CTD{2}.performance.off_axis{2}.resnormS]);
    end
    %----------------------------------------------------------------------
    fprintf(['CTD:\n']);
    fprintf(['Slope in horizontal planes (RMS): ' num2str(mean(sloperms)) ' samples\n']);
    fprintf(['Variance in sagittal planes: ' num2str(mean(sag_var)) ' samples\n']);
    fprintf(['Max deviation in sagittal planes: ' num2str(max(sag_max)) ' samples\n']);
    fprintf(['Resnorm P (on-axis): ' num2str(mean(resnormPon)) ' samples\n']);
    fprintf(['Resnorm P (off-axis):' num2str(mean(resnormPoff)) ' samples\n']);
    fprintf(['Resnorm S (on-axis): ' num2str(mean(resnormSon)) ' samples\n']);
    fprintf(['Resnorm S (off-axis):' num2str(mean(resnormSoff)) ' samples\n']);

    %AGD-------------------------------------------------------------------
    for ii=1:length(results)
        sloperms(ii)=mean(results(ii).AGD{2}.performance.sloperms);
        sag_var(ii)=mean(results(ii).AGD{2}.performance.sag_var);
        sag_max(ii)=max(results(ii).AGD{2}.performance.sag_max);
        resnormPon(ii)=mean([results(ii).AGD{1}.performance.on_axis{1}.resnormP ...
            results(ii).AGD{1}.performance.on_axis{2}.resnormP]);
        resnormPoff(ii)=mean([results(ii).AGD{2}.performance.off_axis{1}.resnormP ...
            results(ii).AGD{2}.performance.off_axis{2}.resnormP]);
        resnormSon(ii)=mean([results(ii).AGD{1}.performance.on_axis{1}.resnormS ...
            results(ii).AGD{1}.performance.on_axis{2}.resnormS]);
        resnormSoff(ii)=mean([results(ii).AGD{2}.performance.off_axis{1}.resnormS ...
            results(ii).AGD{2}.performance.off_axis{2}.resnormS]);
    end
    %----------------------------------------------------------------------
    fprintf(['AGD:\n']);
    fprintf(['Slope in horizontal planes (RMS): ' num2str(mean(sloperms)) ' samples\n']);
    fprintf(['Variance in sagittal planes: ' num2str(mean(sag_var)) ' samples\n']);
    fprintf(['Max deviation in sagittal planes: ' num2str(max(sag_max)) ' samples\n']);
    fprintf(['Resnorm P (on-axis): ' num2str(mean(resnormPon)) ' samples\n']);
    fprintf(['Resnorm P (off-axis):' num2str(mean(resnormPoff)) ' samples\n']);
    fprintf(['Resnorm S (on-axis): ' num2str(mean(resnormSon)) ' samples\n']);
    fprintf(['Resnorm S (off-axis):' num2str(mean(resnormSoff)) ' samples\n']);

    %MCM-------------------------------------------------------------------
    for ii=1:length(results)
        sloperms(ii)=mean(results(ii).MCM{2}.performance.sloperms);
        sag_var(ii)=mean(results(ii).MCM{2}.performance.sag_var);
        sag_max(ii)=max(results(ii).MCM{2}.performance.sag_max);
        resnormPon(ii)=mean([results(ii).MCM{1}.performance.on_axis{1}.resnormP ...
            results(ii).MCM{1}.performance.on_axis{2}.resnormP]);
        resnormPoff(ii)=mean([results(ii).MCM{2}.performance.off_axis{1}.resnormP ...
            results(ii).MCM{2}.performance.off_axis{2}.resnormP]);
        resnormSon(ii)=mean([results(ii).MCM{1}.performance.on_axis{1}.resnormS ...
            results(ii).MCM{1}.performance.on_axis{2}.resnormS]);
        resnormSoff(ii)=mean([results(ii).MCM{2}.performance.off_axis{1}.resnormS ...
            results(ii).MCM{2}.performance.off_axis{2}.resnormS]);
    end
    %--------------------------------------------------------------------------
    fprintf(['MCM:\n']);
    fprintf(['Slope in horizontal planes (RMS): ' num2str(mean(sloperms)) ' samples\n']);
    fprintf(['Variance in sagittal planes: ' num2str(mean(sag_var)) ' samples\n']);
    fprintf(['Max deviation in sagittal planes: ' num2str(max(sag_max)) ' samples\n']);
    fprintf(['Resnorm P (on-axis): ' num2str(mean(resnormPon)) ' samples\n']);
    fprintf(['Resnorm P (off-axis):' num2str(mean(resnormPoff)) ' samples\n']);
    fprintf(['Resnorm S (on-axis): ' num2str(mean(resnormSon)) ' samples\n']);
    fprintf(['Resnorm S (off-axis):' num2str(mean(resnormSoff)) ' samples\n']);
end

end

function idx=ARI_FindPosition(data,azimuth,elevation)
    psi=sin(elevation/180*pi).*sin(data.APV(:,2)/180*pi) + ...
        cos(elevation/180*pi).*cos(data.APV(:,2)/180*pi).*...
        cos(azimuth/180*pi-data.APV(:,1)/180*pi);
    [~,idx]=min(acos(psi));
end