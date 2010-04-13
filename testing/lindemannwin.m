function sig = lindemannwin(sig,N_1)
%LINDEMANNWIN Applies a linear onset window
%   Usage: outsig = lindemannwin(sig,N_1)
%
%   Input parameters:
%       sig - two channel signal to which the window should be apllied
%       N_1 - start time for the binaural cross-correlation (samples)
%
%   Output parameters:
%       sig - the given input signal but with the applied window
%
%   LINDEMANNWIN(sig,N_1) applies a linear onset window of length N_1/2 to the
%   given signal and returns the signal. Lindemann has applied this window to
%   the signals he used in his paper (see lindemann1986a, p. 1614).
%
%   See also: lindemann
%
%R lindemann1986a
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameters ------------------------------------

error(nargchk(2,2,nargin));

if ~isnumeric(sig) || size(sig,2)~=2
    error('%s: sig has to be a signal vector of size n x 2.',upper(mfilename));
end

if ~isnumeric(N_1) || ~isscalar(N_1) || N_1<=0
    error('%s: N_1 has to be a positive scalar.',upper(mfilename));
end


% ------ Computation -----------------------------------------------------
siglen = size(sig,1);
if siglen<ceil(N_1/2)
    error('%s: the length of sig has to be longer than ceil(N_1/2).',...
        upper(mfilename));
end
% Generate window
win = [ linspace(0,1,ceil(N_1/2)) zeros(1,siglen-ceil(N_1/2)) ]';
% Apply window
sig(:,1) = win .* sig(:,1);
sig(:,2) = win .* sig(:,2);
