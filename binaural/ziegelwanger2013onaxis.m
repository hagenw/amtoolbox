function y=ziegelwanger2013onaxis(p,x)
%ZIEGELWANGER2013ONAXIS
%   y=ziegelwanger2013onaxis(p,x)
%
%   Input:
%       p: 
%           r......sphere radius
%           phi....position of the ear (azimuth angle)
%           theta..position of the ear (elevation angle)
%           delay..constant delay
%
%       x:
%
%   Output:
%       y:

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna, Austria



	y=p(1)/343.*( ...
	       (sign(sin(p(3)).*sin(x(:,2))+cos(p(3)).*cos(x(:,2)).*cos(p(2)-x(:,1)))/2+0.5).* ...
	       (1-sin(p(3)).*sin(x(:,2))-cos(p(3)).*cos(x(:,2)).*cos(p(2)-x(:,1)))+ ...
	       (-sign(sin(p(3)).*sin(x(:,2))+cos(p(3)).*cos(x(:,2)).*cos(p(2)-x(:,1)))/2+0.5).* ...
	       (1+acos(sin(p(3)).*sin(x(:,2))+cos(p(3))*cos(x(:,2)).*cos(p(2)-x(:,1)))-pi/2))+p(4)-p(1)/343;
end