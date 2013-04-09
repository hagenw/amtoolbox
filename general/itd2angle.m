function angle = itd2angle(itd,ild,f_inst,tr)
% % % 20.july.12; slight update on the function call for polyfit

load poly_lookup

unwrapped_itd = itd + round(0.4*sign(round(ild/2/(abs(tr)+1e-9)))-0.4*sign(itd))./f_inst;
angle = zeros(size(itd));

for n = 1:size(itd,2)
    %[p ,S ,MU]=polyfit(lookup.mitds(:,n),lookup.azi,9);
    angle(:,n)=polyval(p(:,n),unwrapped_itd(:,n),S{n} ,MU(:,n));
%   by calling the output S and MU, lookup.mitds is z-scored, thus improving the fitting
end
% neglect angles > 90°. WARNING => systematic underestimation for azi ~ 90°
angle(abs(angle)>90) = NaN;
%angle(abs(angle)>90) = 0*sign(angle(abs(angle)>90));
%disp('warning: temporary change')

