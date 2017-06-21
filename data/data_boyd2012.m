function data = data_boyd2012
%DATA_BOYD2012 - Data from Boyd et al. (JASA-EL, 2012)
%   Usage: data = data_boyd2012
%
%   Mean externalization scores of NH listeners extracted from top panels 
%   (1 talker condition) of Fig. 1 
%
%   Output parameters:
%     data    : structure with fields
%                 reference
%                 mix
%                 ITE with subfields BB and LP
%                 BTE with subfields BB and LP
%
%   References: Boyd et al. (JASA-EL, 2012)

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

data.reference = 100;
data.mix = 100:-25:0;
data.ITE.BB = [100,95,58,45,33];
data.BTE.BB = [63,66,38,38,24];
data.ITE.LP = [72,75,42,26,12];
data.BTE.LP = [61,62,33,24,16];

end