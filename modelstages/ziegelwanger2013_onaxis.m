function y=ziegelwanger2013_onaxis(p,x)
%ziegelwanger2013_onaxis XXX
%   Usage: y=ziegelwanger2013_onaxis(p,x)
%
%   Input parameters:
%       p: on-axis model parameters [SI-units]
%       x: HRTF direction (azimuth,elevation) [rad]
%   Output parameters:
%       y: time-of-arrival [s]
%
%   `toa=ziegelwanger2013_onaxis(p,x)` calculates time-of-arrivals (TOAs) for
%   given model parameters (p) and directions (x) with an on-axis
%   time-of-arrival model.
%
%   See also: ziegelwanger2013, ziegelwanger2013_offaxis,
%   data_ziegelwanger2013, exp_ziegelwanger2013
%
%   References: ziegelwanger2013 majdak2013toa

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria
    
r=p(1); %............. sphere radius [m]
phi_ear=p(2); %....... position of the ear (azimuth angle) [rad]
theta_ear=p(3); %..... position of the ear (elevation angle) [rad]
delay=p(4); %......... constant delay [s]

y=r/343.*( ...
       (sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
       (1-sin(theta_ear).*sin(x(:,2))-cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))+ ...
       (-sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
       (1+acos(sin(theta_ear).*sin(x(:,2))+cos(theta_ear)*cos(x(:,2)).*cos(phi_ear-x(:,1)))-pi/2))+delay-r/343;
end