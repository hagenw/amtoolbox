function y=ziegelwanger2014offaxis(p,x)
%ZIEGELWANGER2013OFFAXIS XXX
%   Usage: y=ziegelwanger2014offaxis(p,x) 
%
%   Input:
%       p: off-axis time-of-arrival model parameters [SI-units]
%       x: HRTF direction (azimuth,elevation) [rad]
%   Output:
%       y: time-of-arrival [s]
%
%   `toa=ziegelwanger2014offaxis(p,x)` calculates time-of-arrivals for given
%   model parameters (p) and directions (x) with an off-axis time-of-arrival
%   model.
%
%   See also: ziegelwanger2014, ziegelwanger2014onaxis,
%   data_ziegelwanger2014, exp_ziegelwanger2014
%
%   References: ziegelwanger2014

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria

r=p(1); ............. sphere radius [m]
xM=p(2); ............ x-coordinate of the sphere center [m]
yM=p(3); ............ y-coordinate of the sphere center [m]
zM=p(4); ............ z-coordinate of the sphere center [m]
delay=p(5); ......... constant dely [s]
phi_ear=p(6); ....... position of the ear (azimuth angle) [rad]
theta_ear=p(7); ..... position of the ear (elevation angle) [rad]

M=sqrt(xM^2+yM^2+zM^2);

beta=acos(-cos(x(:,2)).*(xM*cos(x(:,1))+yM*sin(x(:,1)))-zM*sin(x(:,2)));
s2=-r+M*cos(beta)+sqrt(r^2+M^2*cos(beta).^2+2*M*r);
gamma=pi-beta-acos((2*M^2+2*M*r-2*r*s2-s2.^2)/(2*M^2+2*M*r));
if M==0
    s1=zeros(size(x,1),1);
else
    s1=M*cos(beta)./(2*(M+r).*tan(gamma/2));
end

y=1/340*((r* ...
    ((sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
    (1-sin(theta_ear).*sin(x(:,2))-cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))+ ...
    (-sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
    (1+acos(sin(theta_ear).*sin(x(:,2))+cos(theta_ear)*cos(x(:,2)).*cos(phi_ear-x(:,1)))-pi/2))) ...
    +s1+s2) ...
    +delay-(M+r)/340;

end