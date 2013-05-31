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
narginchk(2,2);
if ( ~isnumeric(ii) || ~isscalar(ii) )
    error('%s: ii has to be a scalar.',upper(mfilename));
end
if ( ~isnumeric(nii) || ~isscalar(ii) )
    error('%s: nii has to be a scalar.',upper(mfilename));
end


% ------ Generate the progress bar ---------------------------------------
% \r is not working under some Windows systems, therefore we use \b to clear the
% line
str_length = 30; % this allows to handle integers up to 12 characters
if ii==nii
    clear_line(str_length);
    fprintf(1,'Run %.0f/%.0f\n',ii,nii);
else
    clear_line(str_length);
    fprintf(1,'Run %.0f/%.0f',ii,nii);
end
% Octave didn't show the output directly in a function call, in order to do so
% it has explicitly flushed to stdout
if isoctave
    fflush(1);
end
end % end of function


%% ------ Subfunctions ---------------------------------------------------
function clear_line(str_length)
    for ii=1:str_length
        fprintf(1,'\b');
    end
end
