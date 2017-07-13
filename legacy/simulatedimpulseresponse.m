function  [out,fs]=simulatedimpulseresponse(varargin)

warning('Warning: SIMULATEDIMPULSERESPONSE will be removed in a future release. Use SIG_JOERGENSEN2011 instead. ');  

[out,fs]=sig_joergensen2011(varargin{:});