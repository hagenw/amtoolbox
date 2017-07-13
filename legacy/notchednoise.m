function out=notchednoise(varargin)

warning('Warning: NOTCHEDNOISE will be removed in a future release. Use SIG_NOTCHEDNOISE instead. ');  

out=sig_notchednoise(varargin{:});