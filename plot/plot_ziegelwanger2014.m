function [h]=plot_ziegelwanger2014(Obj,data,type,color,ele,ch,time,style,width)
%plot_ziegelwanger2014   Plot time-of-arrival
%   Usage: plot_ziegelwanger2014(Obj,data,type,color,ele,ch,time,style,width)
%
%   `plot_ziegelwanger2014(Obj,data,type,color,ele,ch,time,style,width)`
%   plots TOA-data in horizontal planes.
%
%   Input:
%       Obj: SOFA object
% 
%       data:
% 
%       type:
%
%       color:
%
%       ele:
%
%       ch:
%
%       time:
%
%       style:
%
%       width:
%
% 
%   Output:
%       h: figure handle
%
%   Estimates the Time-of-Arrival for each column in input data hM and corrects 
%   the results with a geometrical model of the head.
%
%   Examples:
%   ---------
% 
%   To plot the modelled TOA in the horizontal plane after using
%   ziegelwanger2014, use::
%
%       plot_ziegelwanger2014((Obj,Obj.Data.Delay,1,'b',0,1,1);
%
%   See also: ziegelwanger2014_onaxis, ziegelwanger2014_offaxis,
%   data_ziegelwanger2014, exp_ziegelwanger2014
%
%   References: ziegelwanger2014

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria


%% ----------------------------check variables-----------------------------
if time
    data=data/Obj.Data.SamplingRate*1000;
end

if exist('style','var')
    if ischar(style)
        style={style};
    end
end

if ~ischar(color)
    rgbcolor=color;
    color='b'*ones(1,size(color,1));
end

Obj.APV=Obj.SourcePosition(:,1:2);

%% --------------------------------plot TOA--------------------------------
[~,idx1]=sort(Obj.APV(:,1));
idx2=find(Obj.APV(idx1,2)==ele);

for ii=1:length(type)
    switch type(ii)
        case 1
            h=plot(Obj.APV(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) '-']);
        case 2
            h=plot(Obj.APV(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) '--']);
        case 3
            idx3=find(Obj.APV(idx1(idx2),2)==ele & data(idx1(idx2),ii*2-1+(ch-1))~=0);
            if exist('style','var')
                h=plot(Obj.APV(idx1(idx2(idx3)),1),data(idx1(idx2(idx3)),ii*2-1+(ch-1)),[color(ii) style{ii}],'MarkerSize',width(ii));%,'MarkerFaceColor',color(ii)
            else
                h=plot(Obj.APV(idx1(idx2(idx3)),1),data(idx1(idx2(idx3)),ii*2-1+(ch-1)),[color(ii) 'o'],'MarkerSize',2,'MarkerFaceColor',color(ii));
            end
        case 4
            h=plot(Obj.APV(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) style{ii}],'LineWidth',width(ii));
    end
    if exist('rgbcolor','var')
        set(h,'color',rgbcolor(ii,:));
    end
    hold on
end

title(['Elevation: ' num2str(ele) 'deg'])
xlim([0 359])
xlabel('Azimuth (deg)')
if time
    ylabel('Time (ms)')
else
    ylabel('Time (samples)')
end
grid on
set(gca,'XTick',[0 90 180 270 360]);

end %of funciton