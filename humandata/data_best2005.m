function data = data_best2005
%DATA_BEST2005  Listener averages of absolute polar angle error and SCC
%   Usage: data = data_best2005
%
%   Output parameters:
%     data.qe     : quadrant error rates
%     data.ape    : absolute polar angle error
%     data.seape  : standard error of absolute polar angle error
%     data.meta   : condition labels corresponding to data entries
%
%   `data_best2005` returns listener averages of absolute polar angle error
%   and SCCs from Best et al. (2005): Fig.5(b), Fig.10(b), and Tab.2
%
%
%   References best2005highfrequency

% AUTHOR: Robert Baumgartner

% Mean quadrant error rates from Tab.2 (Exp. I)
data.qe =     [3.0,8.4,nan,nan,16.4];

% Mean absolute polar angle errors from Fig.5b (Exp. I) Fig.10b (Exp. II)
data.ape =    [21,30,38,42.5,46];
data.seape =  ones(1,5);

data.meta =   {'BB noise','BB speech','-20dB','-40dB','-60dB'};
end