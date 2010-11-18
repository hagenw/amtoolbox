function [med,pol]=bmp2gr(bmpname)
% BMP2GR converts a bitmap to gain response data (in dB) with conditions
% Usage:  [med,pol]=bmp2gr(bmpname)
% Input arguments:
%       bmpname:    filename of bitmap
% Output arguments:
%       medir:   	gain response data (in dB) on median plane for 
%                   53 polar angles defined in pol
%       pol:     	53 equally spaced polar angles between -55° and 235°
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bmp=imread(bmpname,'bmp');
hrtf=zeros(size(bmp,2),size(bmp,1));

sq=4;
for jj=1:size(hrtf,2)
    if jj<=sq || jj>=size(hrtf,2)-sq
        av=0;
    else
        av=sq;
    end
    for ii=1:size(hrtf,1)
        temp=mean(bmp(end-jj+1-av:end-jj+1+av,ii,1));
        if temp<=25 && temp >=15
            hrtf(ii,jj)=-15;
        elseif temp<=50 && temp>40
            hrtf(ii,jj)=-10;
        elseif temp<=85 && temp>65
            hrtf(ii,jj)=-5;
        elseif temp<=150 && temp>120
            hrtf(ii,jj)=0;
        elseif temp<=185 && temp>160
            hrtf(ii,jj)=5;
        elseif temp<=210 && temp>199
            hrtf(ii,jj)=10;
        elseif temp>220
            hrtf(ii,jj)=15;
        else
            if ii==1
                hrtf(ii,jj)=0;
            else
                hrtf(ii,jj)=hrtf(ii-1,jj);
            end
        end
    end
end

sq=3;
for jj=1:size(hrtf,2)
    if jj<=sq || jj>=size(hrtf,2)-sq
        av=0;
    else
        av=sq;
    end
    for ii=1:size(hrtf,1)
        hrtf(ii,jj)=mean(hrtf(ii,jj-av:jj+av));
    end
end

% figure; pcolor(hrtf'),shading flat
posidx=5:9:475;
pol=-55:290/52:235;
med=hrtf(:,posidx);
% figure; pcolor(med'),shading flat

end