function itd = dietz2011unwrapitd(itd,ild,f_inst,tr)
%DIETZ2011UNWRAPITD unwraps the given itd using the sign of the corresponding ild
%   Usage: itd = dietz2011unwrapitd(itd,ild,f_inst,tr)
%
%   Input parameters:
%       itd    : itd to unwrap
%       ild    : corresponding ild value
%       f_inst : instantaneous frequency
%       tr     : only apply the unwrap mechanism for ild greater than the
%                threshold tr, because for values near 0 it could be wrong
%                (default: 2.5)
%
%   Output parameters:
%       itd    : unwrapped itd
%
%   `dietz2011unwrapitd(itd,ild,f_inst,tr)` unwraps the given itd using the sign of the
%   corresponding ild value. Unwrapping means, that the ild value is used to
%   decide to which direction the itd belongs, which can be unclear for
%   large values, because of the $\pi$ limit (see Dietz et al. 2011, Fig. 2).
%
%   See also: dietz2011, wierstorf2013
%
%   References: dietz2011auditory wierstorf2013

% AUTHOR: Mathias Dietz, Hagen Wierstorf (for AMT)

%% ===== Checking of input parameters ===================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
if nargin==3
    tr = 2.5;
end


%% ===== Calculation ====================================================
itd = itd + ...
    round( ... % this will be -1,0,1
        0.4*sign(round(ild/2 / (abs(tr)+1e-9))) - 0.4*sign(itd) ) ...
    ./ f_inst;
