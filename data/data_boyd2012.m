function data = data_boyd2012
%DATA_BOYD2012 - Data from Boyd et al. (JASA-EL, 2012)
%
%   Usage: data = data_boyd2012
%
%   Mean externalization scores of NH listeners extracted from top panels 
%   (1 talker condition) of Fig. 1 
%
%   Output parameters:
%     data    : structure with fields
%                 ID ... subject ID
%                 Resp ... externalization responses 
%                 BRIR ... binaural room impulse responses for 4 positions
%                 Target ... target stimuli (single and four talker conditions)
%                 Reference_1T ... reference stimuli (single talker)
%                 Reference_4T ... reference stimuli (four talkers)
%                 fs ... sampling rate in Hz
%
%   References: boyd2012

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

tmp = amt_load('boyd2012','boyd2012.mat');
data = tmp.data;

end