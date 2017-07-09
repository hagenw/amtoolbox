function out=bmsin(varargin)

warning('Warning: BMSIN will be removed in a future release. Use SIG_LINDEMANN1986 instead. ');  

out=sig_lindemann1986(varargin{:});