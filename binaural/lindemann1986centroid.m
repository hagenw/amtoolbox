function d = lindemann1986centroid(cc)
%LINDEMANN1986CENTROID Calculates the centroid for a cross-correlation
%   Usage: d = lindemann1986centroid(cc)
%
%   Input parameters:
%       cc  : Lindemann cross-correlation. Dim: 1 x delay line length
%
%   Output parameters:
%       d   : lindemann1986centroid in the range -1..1~ms
%
%   `lindemann1986centroid(cc)` calculates the centroid for a given
%   cross-correlation from the Lindemann model.
%
%   The centroid is computed by (see Lindemann (1986a), page 1613, eq. 22):
%
%   ..         M                  M
%       d = ( sum m*Psi(m) ) / ( sum Psi(m) )
%             m=-M               m=-M
%
%   .. math:: d = \frac{\sum_{m=-M}^{M} m*\Psi (m)}{\sum_{m=-M}^M \Psi (m) }
%
%   where *M* is half the length of the delay line $-M,...,M$.
%
%   See also: lindemann1986
%
%   References: lindemann1986a
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input  parameters -----------------------------------

error(nargchk(1,1,nargin));

if ~isnumeric(cc) || ~isvector(cc)
    error('%s: cc has to be a numeric vector signal!',upper(mfilename));
end
% Ensure size(cc) = delay line length x 1
if size(cc,1)==1
    cc = cc';
end


% ------ Computation -----------------------------------------------------
% Calculate the length of the delay line as -M:M
m = linspace(-1,1,length(cc))';
% Calculate the centroid using the -M:M delay line
d = sum(m.*cc) / sum(cc);

