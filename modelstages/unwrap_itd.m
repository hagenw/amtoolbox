function itd = unwrap_itd(itd,ild,f_inst)
%UNWRAP_ITD unwraps the given itd using the sign of the corresponding ild
%   Usage: itd = unwrap(itd,ild)
%
%   Input parameters:
%       itd    : itd to unwrap
%       ild    : corresponding ild value
%       f_inst : instantaneous frequency
%
%   Output parameters:
%       itd    : unwrapped itd
%
%   `unwrap_itd(itd,ild)` unwraps the given itd using the sign of the
%   corresponding ild value. Unwrapping means, that the ild value is used to
%   decide to which direction the itd belongs, which can be unclear for
%   large values, because of the $\pi$ limit (see Dietz et al. 2011).
%
%   See also: dietz2011
%
%   References: dietz2011auditory

% AUTHOR: Hagen Wierstorf

%% ===== Checking of input parameters ===================================
nargmin = 3;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargnumeric(itd,ild,f_inst);


%% ===== Configuration ==================================================
% Only apply the unwrap mechanism for ild greater than a threshold value,
% because for values near 0 it could be wrong.
ild_threshold = 2.5;


%% ===== Calculation ====================================================
itd = itd + ...
    round( ... % this will be -1,0,1
    0.4*sign(round(ild(:,1:12)/2 / (abs(ild_threshold)+1e-9))) ...
    -0.4*sign(itd)...
    ) ./ f_inst;
