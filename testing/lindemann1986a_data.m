function data = lindemann1986a_data(fignum,flag)
%LINDEMANN1986a_DATA Returns data points from the lindemann papers
%   Usage: data = lindemann1986a_data(fignum,flag)
%          data = lindemann1986a_data(fignum)
%
%   Input parameters:
%       fignum  - number of the figure for which the data should be provided
%       flag    - optinal flag to specify which data from the given figure
%                 should be returned
%
%   Output parameters:
%       data    - the data points from the given figure as a matrix with x and y
%                 data or (in cases where no clear x values were given such as
%                 T/2) as a vector containing only the y data
%
%   LINDEMANN1986a_DATA(fig_num) returns data points for the given figure
%   number from his 1986a paper.
%
%R  lindemann1986a
%
%   See also: testlindemann, lindemanndata1986b
%

%   AUTHOR: Hagen Wierstorf


%% ------ Checking of input  parameters ---------------------------------

error(nargchk(1,2,nargin));

if ~isnumeric(fignum) || ~isscalar(fignum) || fignum<=0
    error('%s: fignum has to be a positive scalar!',upper(mfilename));
end

if nargin == 2 && ~ischar(flag)
    error('%s: flag has to be a string!',upper(mfilename));
end


%% ------ Data points from the paper ------------------------------------
% The following data points are guessed after the plots in the paper from
% Lindemann 1986a.
%
% Data for the given figure
if fignum==11
    % Return data for ILD | lateral displacement (10 is max displacement)
    if nargin~=2
        error(['%s: for fig. 11 you have to specify the data: "yost" or ',...
            '"sayers"!'],upper(mfilename));
    elseif strcmpi(flag,'yost')
        data = [ 0  0.0
                 3  2.1
                 6  4.4
                 9  6.4
                12  8.6
                15  9.2
                18 10.0 ];
    elseif strcmpi(flag,'sayers')
        data = [ 6  3.9
                 9  4.9
                12  8.9 ];
    else
        error('%s: flag has to be "yost" or "sayers" for fig. 11!', ...
            upper(mfilename));
    end
elseif fignum==12
    if nargin~=2
        error(['%s: for fig. 11 you have to specify the data: "400" or ',...
            '"600"!'],upper(mfilename));
    elseif strcmpi(flag,'400')
        data = [ 0.4
                 3.3
                 6.4
                 8.3
                 6.7
                 3.3
                 0.1 ];
    elseif strcmpi(flag,'600')
        data = [ 0.6
                 2.6
                 5.8
                 7.2
                 7.0
                 2.6
                 0.4 ];
    else
        error('%s: flag has to be "400" or "600" for fig. 12!', ...
            upper(mfilename));
    end
elseif fignum==13
    data = [];
else
    error('%s: no data are avaiable for fig. %i!',upper(mfilename),fignum);
end
