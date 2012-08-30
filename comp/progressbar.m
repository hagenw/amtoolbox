function progressbar(ii,nii)
% PROGRESSBAR Show the progress of an iteration
%   Usage: progressbar(ii,nii)
%
%   Input parameters:
%       ii  - Current iteration
%       nii - number of iterations
%
%   PROGRESSBAR(ii,nii) displays the progress of a loop.
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameters ------------------------------------

error(nargchk(2,2,nargin));

if ( ~isnumeric(ii) || ~isscalar(ii) )
    error('%s: ii has to be a scalar.',upper(mfilename));
end

if ( ~isnumeric(nii) || ~isscalar(ii) )
    error('%s: nii has to be a scalar.',upper(mfilename));
end


% ------ Generate the progress bar ---------------------------------------

if ii==nii
    fprintf(1,'\rRun %.0f/%.0f\n',ii,nii);
else
    fprintf(1,'\rRun %.0f/%.0f',ii,nii);
end

% Octave didn't show the output directly in a function call, in order to do so
% it has explicitly flushed to stdout
if isoctave
    fflush(1);
end
