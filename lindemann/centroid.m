function d = centroid(cc)
%CENTROID Calculates the centroid for a cross-correlation
%   Usage: d = centroid(cc)
%
%   Input parameters:
%       cc  - Lindemann cross-correlation. Dim: 1 x delay line length
%
%   Output parameters:
%       d   - centroid in the range -1..1~ms
%
%   CENTROID(cc) calculates the centroid for a given cross-correlation from the
%   Lindemann model.
%
%   The centroid is computed by (see lindemann1986a, page 1613, eq. 22):
%
%C             M                  M
%C      d = ( sum m*Psi(m) ) / ( sum Psi(m) )
%C            m=-M               m=-M
%
%   where M is half the length of the delay line -M:M.
%
%   See also: lindemann
%
%R lindemann1986a
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input  parameters -----------------------------------

error(nargchk(1,1,nargin));

if ~isnumeric(cc) || ~isvector(cc)
    error('%s: cc has to be a numeric vector signal!',upper(mfilename));
end


% ------ Computation -----------------------------------------------------
% Calculate the length of the delay line as -M:M
M = fix(length(cc)/2);
% Calculate the centroid using the -M:M delay line
d = sum((-M:M)'/M.*cc)/sum(cc);
