function data = data_hassager2016
%DATA_HASSAGER2016 - Data from Hassager et al. (JASA, 2016)
%   Usage: data = data_hassager2016
%
%   NH data from Hassager et al. (JASA, 2016), Fig. 6,
%   representing listeners' ratings for the dir(ect-sound)
%   condition
%
%   Output parameters:
%     data    : structure with fields
%                 B ... Bandwidth Factor (ERB)
%                 angle ... source angle (deg)
%                 rating ... Externalization Rating (dim: B x angle)
%
%   References: Hassager et al. (JASA, 2016)

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

data.B = [nan,0.316, 0.570, 1.03, 1.85, 3.33, 6.0, 10.8, 19.5 35.0, 63.1];

data.angle = [0,50];

data.rating = ...
  [ 4.79,4.73,4.70,4.67,4.43,3.65,2.81,1.97,1.60,1.49,1.30;
    4.94,4.92,4.94,4.85,4.79,4.33,3.86,3.21,2.59,2.08,1.60 ]';

end