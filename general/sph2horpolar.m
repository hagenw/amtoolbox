function [lat,pol]=sph2horpolar(azi,ele)
% SPH2HORPOLAR from spherical to horizontal-polar coordinate system
%
% Usage:    [lat,pol] = sph2horpolar(azi,ele)
%
% Input parameters:
%     azi     : azimuth in deg
%     ele     : elevation in deg
%
% Output parameters:
%     lat     : lateral angle in deg, [-90°,+90°]
%     pol     : polar angle in deg, [-90°,270°]
%
%   `geo2horpolar(...)` converts spherical coordinates (azimuth and 
%   elevation) into the coordinates of the horizontal-polar system, i.e.,
%   lateral and polar angles.
%
% AUTHOR: Piotr Majdak, Robert Baumgartner

azi=mod(azi+360,360);
ele=mod(ele+360,360);

razi = deg2rad(azi);
rele = deg2rad(ele);
rlat=asin(sin(razi).*cos(rele));
rpol=asin(sin(rele)./cos(rlat));
rpol(cos(rlat)==0)=0;
pol = rad2deg(rpol);
lat = rad2deg(rlat);

idx = find(razi>pi/2 & razi < 3*pi/2 & (rele < pi/2 | rele > 3*pi/2));
pol(idx)=180-pol(idx);
idx = find(~(razi>pi/2 & razi < 3*pi/2) & rele > pi/2 & rele < 3*pi/2);
pol(idx)=180-pol(idx);

end