function data = data_middlebrooks1999
%DATA_MIDDLEBROOKS1999 Statistics about non-individualized HRTFs
%   Usage: data = data_middlebrooks1999
%
%   Output parameters:
%     data.qe_own   : quadrant error rate (QE) when localizing with own
%                     HRTFs
%     data.qe_other : QE when localizing with others' HRTFs
%     data.pe_own   : local polar RMS error (PE) when localizing with own
%                     HRTFs
%     data.pe_other : PE when localizing with others' HRTFs
%     data.pb_own   : magnitude of polar bias (PB) when localizing with own
%                     HRTFs; upper-rear quadrant excluded from analysis
%     data.pb_other : PB when localizing with others' HRTFs
%
%   `data_middlebrooks1999` returns statistics summary from Fig. 13
%   (Middlebrooks, 1999b) showing the effect of non-individualized HRTFs.
%
%   Statistics of those parameters are stored as *.mean* and *.quantiles*
%   representing the arithmetic mean and {0,5,25,50,75,95,100} quantiles,
%   respectively.
%
%   References: middlebrooks1999nonindividualized

% AUTHOR: Robert Baumgartner

% Quantiles: {0,5,25,50,75,95,100}%
  % QE
  data.qe_own.quantiles = [0,0,1,3,5,13,17];
  data.qe_own.mean = 3;
  data.qe_other.quantiles = [7.5,7.5,13,19,28,38,39];
  data.qe_other.mean = 21;
  % PE
  data.pe_own.quantiles = [22,23,25,27,30,34,36];
  data.pe_own.mean = 28;
  data.pe_other.quantiles = [23,33,38,42,48,54,55];
  data.pe_other.mean = 43;
  % EB
  data.pb_own.quantiles = [1,2.5,6,10,13,20,25.5];
  data.pb_own.mean = 10;
  data.pb_other.quantiles = [0.5,2,6.5,18,29,42,52];
  data.pb_other.mean = 19;
end