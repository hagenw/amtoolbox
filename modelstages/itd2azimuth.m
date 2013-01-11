function phi = itd2azimuth(itd,lookup)
%ITD2AZIMUTH estimates the azimuth for the given itd
%
%   Usage: phi = itd2azimuth(itd,lookup);
%
%   Input parameters:
%       itd    - interaural time difference, ITD (s)
%       lookup - lookup table to use
%
%   Output parameters:
%       phi    - azimuth angle (rad)
%
%   ITD2AZIMUTH(ITD,LOOKUP) estimates the azimuth angle by the given
%   interaural time difference (itd) and the given lookup table.
%   Note: only the first 12 auditory channels are considered.
%
%   see also: lookup_table, ild2azimuth
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ===================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(itd);
isargstruct(lookup);


%% ===== Computation ====================================================
% Fit the lookup data
p = zeros(12,13);
for n = 1:12
    p(n,:)=polyfit(lookup.itd(:,n),lookup.azimuth',12);
end
% Use the fit parameter to estimate the azimuth angle phi with the
% calculated ITD
phi = zeros(size(itd,1),12);
for n = 1:12
    phi(:,n) = polyval(p(n,:),itd(:,n));
end
% Check if we got azimuth values, that extend our range
phi(abs(phi)>rad(95)) = NaN;
