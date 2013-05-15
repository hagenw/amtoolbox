function phi = itd2angle(itd,lookup)
% ITD2ANGLE converts the given ITD to an angle using a lookup table
%   Usage: phi = itd2angle(itd,ild,f_inst,tr)
%
%   Input parameters:
%       itd     : ITDs to convert to angles
%       lookup  : a struct containing the polinomial fitting entries p,MU,S.
%                 This struct can be generated by `itd2anglelookuptable`
%
%   Output parameters:
%       phi     : angles for the corresponding ITD values / deg
%
%   `itd2angle(itd,lookup)` converts the given ITD values to azimuth angles phi.
%   Therefore a lookup table containing the polynomial fitting parameters is
%   used. This lookup table is created from a set of HRTFs with known azimuth
%   angles and stores the corresponding ITD values. Note that itd2angle works
%   with the dietz2011 and with lindemann1986.

% AUTHOR: Mathias Dietz, Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Computation =====================================================
phi = zeros(size(itd));

for n = 1:size(itd,2)
    % by calling the output S and MU, phi is z-scored, thus improving the fitting
    phi(:,n)=polyval(lookup.p(:,n),itd(:,n),lookup.S{n},lookup.MU(:,n));
end
% neglect angles > 95°. WARNING => maybe systematic underestimation for azi ~ 90°
phi(abs(angle)>95) = NaN;

