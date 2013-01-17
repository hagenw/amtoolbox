function [h]=plotziegelwanger2013(data,type,color,ele,ch,time,meta,stimPar,style,width)
%PLOTZIEGELWANGER2013
%   plotziegelwanger2013(data,type,color,ele,ch,time,meta,stimPar,style,width)
%
%   Input:
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
%       meta:
%
%       stimPar:
%
%       style:
%
%       width:
%
%   Output:

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna, Austria


%% ----------------------------check variables-----------------------------
if time
    data=data/stimPar.SamplingRate*1000;
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

%% --------------------------------plot TOA--------------------------------
[~,idx1]=sort(meta.pos(:,1));
idx2=find(meta.pos(idx1,2)==ele);

for ii=1:length(type)
    switch type(ii)
        case 1
            h=plot(meta.pos(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) '-']);
        case 2
            h=plot(meta.pos(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) '--']);
        case 3
            idx3=find(meta.pos(idx1(idx2),2)==ele & data(idx1(idx2),ii*2-1+(ch-1))~=0);
            if exist('style','var')
                h=plot(meta.pos(idx1(idx2(idx3)),1),data(idx1(idx2(idx3)),ii*2-1+(ch-1)),[color(ii) style{ii}],'MarkerSize',width(ii));%,'MarkerFaceColor',color(ii)
            else
                h=plot(meta.pos(idx1(idx2(idx3)),1),data(idx1(idx2(idx3)),ii*2-1+(ch-1)),[color(ii) '.']);
            end
        case 4
            h=plot(meta.pos(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) style{ii}],'LineWidth',width(ii));
    end
    if exist('rgbcolor','var')
        set(h,'color',rgbcolor(ii,:));
    end
    hold on
end

title(['Elevation: ' num2str(ele) '°'])
xlim([0 359])
xlabel('Azimuth (°)')
if time
    ylabel('Time (ms)')
else
    ylabel('Time (samples)')
end
grid on
set(gca,'XTick',[0 90 180 270 360]);

end %of funciton