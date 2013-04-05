function y=ziegelwanger2013onaxis(p,x)
%ZIEGELWANGER2013ONAXIS
%   y=ziegelwanger2013onaxis(p,x) calculates time-of-arrivals (TOAs) for
%   given model parameters (p) and directions (x) with an on-axis
%   time-of-arrival model.
%
%   Input:
%       p: on-axis model parameters [SI-units]
%       x: HRTF direction (azimuth,elevation) [rad]
%   Output:
%       y: time-of-arrival [s]
% 
%   Examples:
%   ---------
% 
%   To calculate TOAs for given parameters p and directions x, use:::
%
%     toa=ziegelwanger2013onaxis(p,x);
%
%   See also: ziegelwanger2013, ziegelwanger2013offaxis,
%   data_ziegelwanger2013, exp_ziegelwanger2013
%
%   References: ziegelwanger2013

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria
    
r=p(1); ............. sphere radius [m]
phi_ear=p(2); ....... position of the ear (azimuth angle) [rad]
theta_ear=p(3); ..... position of the ear (elevation angle) [rad]
delay=p(4); ......... constant delay [s]

y=r/343.*( ...
       (sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
       (1-sin(theta_ear).*sin(x(:,2))-cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))+ ...
       (-sign(sin(theta_ear).*sin(x(:,2))+cos(theta_ear).*cos(x(:,2)).*cos(phi_ear-x(:,1)))/2+0.5).* ...
       (1+acos(sin(theta_ear).*sin(x(:,2))+cos(theta_ear)*cos(x(:,2)).*cos(phi_ear-x(:,1)))-pi/2))+delay-r/343;
end