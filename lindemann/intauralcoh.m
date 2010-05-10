function coher = intauralcoh(insig1,insig2)
% INTAURALCOH Computes the interaural coherence
%   Usage: coherence = intauralcoh(insig1,insig2)
%
%   Input parameters:
%       insig1  - signal vector 
%       insig2  - signal vector of same size as insig1
%
%   Output parameters:
%       coher   - coherence between the two signals
%
%   INTAURALCOH(insig1,insig2) computes the interaural coherence between the two
%   input signals insig1 and insig2 (left and right signal). See lindemann1986a,
%   eq. 26.
%
%R lindemann1986a

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameters ------------------------------------

error(nargchk(2,2,nargin));

if ( ~isnumeric(insig1) || ~isvector(insig1) )
    error('%s: insig1 has to be a vector signal.',upper(mfilename));
end

if ( ~isnumeric(insig2) || ~isvector(insig2) )
    error('%s: insig2 has to be a vector signal.',upper(mfilename));
end


% ------ Computation -----------------------------------------------------

% Check the length of the two input signals
nsig1 = length(insig1);
nsig2 = length(insig2);
if nsig1~=nsig2
    error('%s: insig1 and insig2 has to be of the same length.',...
        upper(mfilename));
end

% Calculate the coherence after lindemann1986a, eq. 26
coherence = zeros(nsig1,1);
div = sqrt(sum(insig1.^2) * sum(insig2.^2));
for tau=1:nsig1
    coherence(tau) = max( sum(insig1.*[insig2(tau:end); insig2(1:tau-1)]) / div );
end
coher = max(coherence);

