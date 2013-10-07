function [Obj,results]=ziegelwanger2014(Obj,estimation,outlierDetection,model,p0_onaxis)
%ZIEGELWANGER2013 Time of arrival estimates
%   usage: meta=ziegelwanger2014(data,estimation,correct,p0_onaxis) 
%
%   Estimates the Time-of-Arrival for each measurement in Obj (SOFA) and
%   corrects the results with a geometrical model of the head.
%
%   Input:
%       Obj: SOFA object
% 
%       estimation (optional): select one of the estimation methods
%           1: Threshold-Detection
%           2: Centroid of squared IR
%           3: Mean Groupdelay
%           4: Minimal-Phase Cross-Correlation (Max) (default)
%           5: Minimal-Phase Cross-Correlation (Centroid)
%           6: Zero-Crossing
%
%       model (optional): correct estimated toa, using geometrical TOA-Model
%           0: TOA estimated
%           1: TOA modeled (default)
%
%       p0_onaxis (optional): startvalues for lsqcurvefit
%           dim 1: [sphere-radius in m,
%                 azimut of ear in radiants,
%                 elevation of ear in radiants, 
%                 direction-independent delay in seconds]
%           dim 2: each record channel
% 
%   Output:
%       Obj: SOFA Object
% 
%       results.toa: data matrix with time of arrival (TOA) for each impulse response (IR):
%           dim 1: each toa in samples
%           dim 2: each record channel
%       results.p_onaxis: estimated on-axis model-parameters
%           dim 1: [sphere-radius in m,
%                 azimut of ear in radiants,
%                 elevation of ear in radiants,
%                 direction-independent delay in seconds]
%           dim 2: each record channel
%       results.p_offaxis: estimated off-axis model-parameters
%           dim 1: [sphere-radius in m,
%                 xM in m,
%                 yM in m,
%                 zM in m,
%                 direction-independent delay in seconds,
%                 azimut of ear in radiants,
%                 elevation of ear in radiants]
%           dim 2: each record channel
%
%   Examples:
%   ---------
% 
%   To calculate the model parameters for the on-axis time-of-arrival model
%   (p_onaxis) and for the off-axis time-of-arrival model (p_offaxis) for a
%   given HRTF set (SOFA object, 'Obj') with the minimum-phase
%   cross-correlation estimation, use::
%
%       [Obj,results]=ziegelwanger2014(Obj,4,1);
%
%   See also: ziegelwanger2014onaxis, ziegelwanger2014offaxis,
%   data_ziegelwanger2014, exp_ziegelwanger2014
%
%   References: ziegelwanger2014

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria

%% ----------------------------convert to SOFA-----------------------------
if ~isfield(Obj,'GLOBAL_Version')
    Obj=SOFAconvertARI2SOFA(Obj.hM,Obj.meta,Obj.stimPar);
end

%% ----------------------------check variables-----------------------------
if ~exist('estimation','var')
    estimation=4;
else if isempty(estimation)
        estimation=4;
    end
end

if ~exist('outlierDetection','var')
    outlierDetection=1;
else if isempty(outlierDetection)
        outlierDetection=1;
    end
end

if ~exist('model','var')
    model=1;
else if isempty(model)
        model=1;
    end
end

if ~exist('p0_onaxis','var')
    p0_onaxis=[[0.0875; pi/2; 0; 0] [0.0875; -pi/2; 0; 0]];
end

%% -------------------------initialize variables---------------------------
p0_onaxis=transpose(p0_onaxis);
p_onaxis=zeros(size(p0_onaxis));
p0_offaxis=zeros(2,7);
p_offaxis=p0_offaxis;

toa=zeros(Obj.API.M,Obj.API.R);
toaEst=zeros(Obj.API.M,Obj.API.R);
indicator=zeros(Obj.API.M,Obj.API.R);
indicator_hor=indicator;
indicator_sag=indicator;

pos(:,1:2)=Obj.SourcePosition(:,1:2);
[pos(:,6),pos(:,7)]=sph2hor(Obj.SourcePosition(:,1),Obj.SourcePosition(:,2));
pos(:,8)=cumsum(ones(Obj.API.M,1));

rs=1;
[hM,hMmin]=ARI_MinimalPhase(Obj,rs);

%% -----------------------estimate time-of-arrival-------------------------
switch estimation
    case 1 %---------------------------Threshold---------------------------
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                toaEst(ii,jj)=find(abs(hM(ii,jj,:))==max(abs(hM(ii,jj,:))),1);
            end
        end
    case 2 %---------------------------Centroid----------------------------
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                toaEst(ii,jj)=find(cumsum(hM(ii,jj,:).^2)>(sum(hM(ii,jj,:).^2)/2),1);
            end
        end
    case 3 %---------------------------Groupdelay--------------------------
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                [Gd,F]=grpdelay(transpose(double(squeeze(hM(ii,jj,:)))),1,Obj.API.N*4,Obj.Data.SamplingRate*4);
                toaEst(ii,jj)=median(Gd(find(F>500):find(F>2000)));
            end
        end
    case 4 %---------------------------Minimal-Phase-----------------------
        corrcoeff=zeros(Obj.API.M,Obj.API.R);
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                [c,lag]=xcorr(squeeze(hM(ii,jj,:)),squeeze(hMmin(ii,jj,:)),Obj.API.N*4-1,'none');
                [corrcoeff(ii,jj),idx]=max(abs(c));
                corrcoeff(ii,jj)=corrcoeff(ii,jj)/sum(hM(ii,jj,:).^2);
                toaEst(ii,jj)=lag(idx);
            end
        end
end
toaEst=toaEst/rs;

%% --------------------Detect-Outliers-in-estimated-TOA--------------------
if outlierDetection
    for ch=1:Obj.API.R
        if outlierDetection==1 || outlierDetection==2
            % Outlier detection: smooth TOA in horizontal planes
            [~,idxSortHor]=sort(pos(:,1));
            epsilon=5;
            slope=zeros(Obj.API.M,1);
            sloperms=[];
            for ele=min(pos(:,2)):epsilon:max(pos(:,2)) %calculate slope for each elevation along azimuth
                idx=find(pos(idxSortHor,2)>ele-epsilon/2 & pos(idxSortHor,2)<=ele+epsilon/2);
                if numel(idx)>1
                    idx(length(idx)+1)=idx(1);
                    slope(idxSortHor(idx(1:end-1)),1)=diff(toaEst(idxSortHor(idx),ch))./(abs(abs(abs(diff(pos(idxSortHor(idx),1)))-180)-180)+0.1);
                    sloperms(end+1)=sqrt(sum(slope(idxSortHor(idx)).^2)/length(slope(idxSortHor(idx))));
                end
            end
%             idx=find(pos(:,2)==0);
%             sloperms=sqrt(sum(slope(idx).^2)/length(slope(idx)))*1.4;
            sloperms=max(sloperms);
            performance.sloperms(ch)=sloperms;
            performance.slopemax(ch)=max(abs(slope(:)));
            for ele=min(pos(:,2)):epsilon:max(pos(:,2))
                idx=find(pos(idxSortHor,2)>ele-epsilon/2 & pos(idxSortHor,2)<=ele+epsilon/2);
                for ii=1:length(idx)
                    if abs(slope(idxSortHor(idx(ii))))>sloperms*2
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
            indicator(:,ch)=indicator_hor(:,ch);
        end
        
        if outlierDetection==1 || outlierDetection==3
            % Outlier detection: constant TOA in sagittal planes
            epsilon=2;
            sag_dev=zeros(Obj.API.M,1);
    %         sag_std=sag_dev;
            for lat=-90:epsilon:90
                idx=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2);
                if length(idx)>2
                    sag_dev(idx,1)=toaEst(idx,ch)-mean(toaEst(idx,ch));
    %                 sag_std(idx,1)=sqrt(sum(sag_dev(idx,1).^2)/length(sag_dev(idx,1)))*ones(length(idx),1);
                end
            end
            sag_std=sqrt(sum(sag_dev.^2)/length(sag_dev));
%             if sag_std<0.6
%                 sag_std=0.6;
%             end
            performance.sag_max(ch)=max(abs(sag_dev));
            performance.sag_std(ch)=sag_std;
            indicator_sag(:,ch)=abs(sag_dev)>2*sag_std;
            indicator(:,ch)=(abs(sag_dev)>2*sag_std | indicator_hor(:,ch));
            clear sag_dev;
        end
    end
end

%% ----------------------Fit-Models-to-estimated-TOA-----------------------
if model
    % Fit on-axis model to outlier adjusted set of estimated TOAs
    for ch=1:Obj.API.R
        p0_onaxis(ch,4)=min(toaEst(indicator(:,ch)==0,ch))/Obj.Data.SamplingRate;
        p0offset_onaxis=[0.06 pi pi 0.001];

        idx=find(indicator(:,ch)==0);
        x=pos(idx,1:2)*pi/180;
        y=toaEst(idx,ch)/Obj.Data.SamplingRate;
        if isoctave
            [~,p_onaxis(ch,:)]=leasqr(x,y,p0_onaxis(ch,:),@ziegelwanger2014onaxis);
        else
            [p_onaxis(ch,:),performance.on_axis{ch}.resnormS,performance.on_axis{ch}.residualS,performance.on_axis{ch}.exitflag,performance.on_axis{ch}.output]=...
                lsqcurvefit(@ziegelwanger2014onaxis,p0_onaxis(ch,:),x,y,p0_onaxis(ch,:)-p0offset_onaxis,p0_onaxis(ch,:)+p0offset_onaxis,optimset('Display','off','TolFun',1e-6));
        end
        toa(:,ch)=ziegelwanger2014onaxis(p_onaxis(ch,:),pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;
    end
    performance.on_axis{1}.resnormS=performance.on_axis{1}.resnormS*Obj.Data.SamplingRate^2/(Obj.API.M-sum(indicator(:,1)));
    performance.on_axis{2}.resnormS=performance.on_axis{2}.resnormS*Obj.Data.SamplingRate^2/(Obj.API.M-sum(indicator(:,1)));
    performance.on_axis{1}.resnormP=norm((toaEst(:,1)-toa(:,1))/Obj.Data.SamplingRate,2)^2*Obj.Data.SamplingRate^2/Obj.API.M;
    performance.on_axis{2}.resnormP=norm((toaEst(:,2)-toa(:,2))/Obj.Data.SamplingRate,2)^2*Obj.Data.SamplingRate^2/Obj.API.M;
    performance.on_axis{1}.residualS=performance.on_axis{1}.residualS*Obj.Data.SamplingRate;
    performance.on_axis{2}.residualS=performance.on_axis{2}.residualS*Obj.Data.SamplingRate;

    % Fit off-axis model to outlier adjusted set of estimated TOAs
    p0_offaxis(1,:)=[mean(p0_onaxis(:,1)) 0 0 0 mean(p0_onaxis(:,4)) p0_onaxis(1,2) p0_onaxis(1,3)];
    p0_offaxis(2,:)=[mean(p0_onaxis(:,1)) 0 0 0 mean(p0_onaxis(:,4)) p0_onaxis(2,2) p0_onaxis(2,3)];
    TolFun=[1e-5; 5e-7];
    for ii=1:size(TolFun,1)
        for ch=1:Obj.API.R
            idx=find(indicator(:,ch)==0);
            x=pos(idx,1:2)*pi/180;
            y=toaEst(idx,ch)/Obj.Data.SamplingRate;
            p0_offaxis(ch,:)=[mean(p0_offaxis(:,1)) mean(p0_offaxis(:,2)) mean(p0_offaxis(:,3)) mean(p0_offaxis(:,4)) mean(p0_offaxis(:,5)) p0_offaxis(ch,6) p0_offaxis(ch,7)];
            p0offset_offaxis=[0.07/ii 0.07/ii 0.07/ii 0.07/ii 0.001 pi pi];
            if isoctave
                [~,p_offaxis(ch,:)]=leasqr(x,y,p0_offaxis(ch,:),@ziegelwanger2014offaxis);
            else
                [p_offaxis(ch,:),performance.off_axis{ch}.resnormS,performance.off_axis{ch}.residualS,performance.off_axis{ch}.exitflag,performance.off_axis{ch}.output]=...
                    lsqcurvefit(@ziegelwanger2014offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0offset_offaxis,p0_offaxis(ch,:)+p0offset_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));
            end
        end
        tmp_p_offaxis=p_offaxis;
        if abs(diff(p_offaxis(:,1)))>0.002 || abs(diff(p_offaxis(:,3)))>0.002
            p_offaxis(:,[1 3])=p_offaxis([2 1],[1 3]);
            for ch=1:Obj.API.R
                idx=find(indicator(:,ch)==0);
                x=pos(idx,1:2)*pi/180;
                y=toaEst(idx,ch)/Obj.Data.SamplingRate;
                p0_offaxis(ch,:)=[p_offaxis(ch,1) mean(p_offaxis(:,2)) p_offaxis(ch,3) mean(p_offaxis(:,4)) mean(p_offaxis(:,5)) p_offaxis(ch,6) p_offaxis(ch,7)];
                p0offset_offaxis=[0.07/ii 0.07/ii 0.07/ii 0.07/ii 0.001 pi/2 pi/2];
                if isoctave
                    [~,p_offaxis(ch,:)]=leasqr(x,y,p0_offaxis(ch,:),@ziegelwanger2014offaxis);
                else
                    [p_offaxis(ch,:),performance.off_axis{ch}.resnormS,performance.off_axis{ch}.residualS,performance.off_axis{ch}.exitflag,performance.off_axis{ch}.output]=...
                        lsqcurvefit(@ziegelwanger2014offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0offset_offaxis,p0_offaxis(ch,:)+p0offset_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));
                end
            end
        end
        if abs(diff(p_offaxis(:,1)))>abs(diff(tmp_p_offaxis(:,1))) && abs(diff(p_offaxis(:,3)))>abs(diff(tmp_p_offaxis(:,3)))
            p_offaxis=tmp_p_offaxis;
        end
        if abs(diff(p_offaxis(:,1)))<0.002 && abs(diff(p_offaxis(:,2)))<0.002 && abs(diff(p_offaxis(:,3)))<0.002 && abs(diff(p_offaxis(:,4)))<0.002
            break
        end
    end
    toa(:,ch)=ziegelwanger2014offaxis(p_offaxis(ch,:),pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;
    performance.off_axis{1}.resnormS=performance.off_axis{1}.resnormS*Obj.Data.SamplingRate^2/(Obj.API.M-sum(indicator(:,1)));
    performance.off_axis{2}.resnormS=performance.off_axis{2}.resnormS*Obj.Data.SamplingRate^2/(Obj.API.M-sum(indicator(:,2)));
    performance.off_axis{1}.resnormP=norm((toaEst(:,1)-toa(:,1))/Obj.Data.SamplingRate,2)^2*Obj.Data.SamplingRate^2/Obj.API.M;
    performance.off_axis{2}.resnormP=norm((toaEst(:,2)-toa(:,2))/Obj.Data.SamplingRate,2)^2*Obj.Data.SamplingRate^2/Obj.API.M;
    performance.off_axis{1}.residualS=performance.off_axis{1}.residualS*Obj.Data.SamplingRate;
    performance.off_axis{2}.residualS=performance.off_axis{2}.residualS*Obj.Data.SamplingRate;
else
    toa=toaEst;
    p_offaxis=p0_offaxis;
end

%Save to output variables
Obj.Data.Delay=toa;
performance.outliers=indicator;
performance.outlierRate=sum(sum(indicator))/Obj.API.M/2*100;
performance.outlierRateL=sum(indicator(:,1))/Obj.API.M*100;
performance.outlierRateR=sum(indicator(:,2))/Obj.API.M*100;
results.toa=toa;
results.p_onaxis=transpose(p_onaxis);
results.p_offaxis=transpose(p_offaxis);
results.performance=performance;
if exist('corrcoeff','var')
    results.performance.corrcoeff=corrcoeff;
end

% for ch=1:Obj.API.R
%     ii=0;
%     count=0;
%     for ele=transpose(unique(Obj.SourcePosition(:,2)))
%         ii=ii+1;
%         if ii==9
%             ii=1;
%         end
%         if ii==1
%             figure('Position',[20 20 600 900]);
%             count=count+1;
%         end
%         subplot(4,2,ii);
%         plotziegelwanger2014(Obj,toaEst,2,'k',ele,ch,0);
%         plotziegelwanger2014(Obj,indicator.*toaEst,3,'r',ele,ch,0);
%         plotziegelwanger2014(Obj,toa,1,'g',ele,ch,0);
% %             ylim([min(toaLIN(:,ch)-sync)/stimPar.SamplingRate*1000-0.5 max(toaLIN(:,ch)-sync)/stimPar.SamplingRate*1000+0.1])
%         xlabel('')
%         ylabel('')
%         if count<3
%             if ii==7 || ii==8
%                 xlabel('Azimut in Grad')
%             end
%         else
%             if ii==5 || ii==6
%                 xlabel('Azimut in Grad')
%             end
%         end
%         if ii==1 || ii==3 || ii==5 || ii==7
%             ylabel('Zeit in ms')
%         end
%     end
% end

end %of function

function [hM,hMmin]=ARI_MinimalPhase(Obj,rs)
    hM=zeros(size(Obj.Data.IR,1),size(Obj.Data.IR,2),size(Obj.Data.IR,3)*rs);
    hMmin=hM;
    for ii=1:Obj.API.M
        for jj=1:Obj.API.R
            hM(ii,jj,:)=resample(Obj.Data.IR(ii,jj,:),rs,1);
        end
    end

    for jj=1:Obj.API.R
        for ii=1:Obj.API.M
%             h=squeeze(Obj.Data.IR(ii,jj,:));
            h=[squeeze(hM(ii,jj,:)); zeros(4096*4-size(hM,3),1)];
            % decompose signal
            amp1=abs(fft(h));

            % transform
            amp2=amp1;
            an2u=-imag(hilbert(log(amp1))); % minimal phase

            % reconstruct signal from amp2 and an2u
            % build a symmetrical phase 
            an2u=an2u(1:floor(length(h)/2)+1);
            an3u=[an2u; -flipud(an2u(2:end+mod(length(h),2)-1))];
            an3=an3u-round(an3u/2/pi)*2*pi;  % wrap around +/-pi: wrap(x)=x-round(x/2/pi)*2*pi
            % amplitude
            amp2=amp2(1:floor(length(h)/2)+1);
            amp3=[amp2; flipud(amp2(2:end+mod(length(h),2)-1))];
            % back to time domain
            h2=real(ifft(amp3.*exp(1i*an3)));
            hMmin(ii,jj,:)=h2(1:Obj.API.N*rs);
        end
    end
end