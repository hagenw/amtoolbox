function angle = itd2angle(itd,ild,f_inst,lookup,tr)

% tr: threshold in dB for assuming no unwrapping (2-3 dB !!!)

p = zeros(size(itd,2),10);
for n = 1:size(itd,2)
    p(n,:)=polyfit(lookup.mitds(:,n),lookup.azi,9);
end

unwrapped_itd = itd + round(0.4*sign(round(ild/2/(abs(tr)+1e-9)))-0.4*sign(itd))./f_inst;
angle = zeros(size(itd));
for n = 1:size(itd,2)
    angle(:,n)=polyval(p(n,:),unwrapped_itd(:,n));
end
% neglect angles > 90°. WARNING => systematic underestimation for azi ~ 90°
angle(abs(angle)>90) = NaN;
%angle(abs(angle)>90) = 0*sign(angle(abs(angle)>90));
%disp('warning: temporary change')

%OLDFORMAT
