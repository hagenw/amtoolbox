function outsig = corpinknoise(siglen,coher)
% CORPINKNOISE Generates pink noise with interaural correlation
%   Usage: outsig = corpinknoise(siglen,coher)
%
%   Input parameters:
%       siglen    - Number of samples of outsig
%       coher       - Interaural coherence of the produced outsig
%
%   Output parameters:
%       outsig  - nsig x 2 correlated pink noise signal

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameters ------------------------------------

error(nargchk(2,2,nargin));

if ( ~isnumeric(siglen) || ~isscalar(siglen) || siglen<0 )
    error('%s: siglen has to be a positive scalar.',upper(mfilename));
end

if ( ~isnumeric(coher) || ~isscalar(coher) || coher<0)
    error('%s: coher has to be a positive scalar.',upper(mfilename));
end


% ------ Computation -----------------------------------------------------

% Generate correlation matrix
R = [1 coher; coher 1];
% Eigen decomposition
[V,D] = eig(R);

% Form correlating filter
W = V*sqrt(D);

% Generate pink noise signal
n = noise(siglen,2,'pink');

% Generate correlated noise
outsig = n * W';
