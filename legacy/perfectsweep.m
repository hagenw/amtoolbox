function out=perfectsweep(varargin)

warning('Warning: PERFECTSWEEP will be removed in a future release. Use SIG_LINSWEEP instead. ');  

out=sig_linsweep(varargin{:});