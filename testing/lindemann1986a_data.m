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
% --- FIG 11 ---
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
% --- FIG 12 ---
elseif fignum==12
    if nargin~=2
        error(['%s: for fig. 11 you have to specify the data: "400" or ',...
            '"600"!'],upper(mfilename));
    elseif strcmpi(flag,'400')
        data = [ -1.00  0.4
                 -0.83  3.3
                 -0.67  6.4
                 -0.50  8.3
                 -0.33  6.7
                 -0.17  3.3
                  0.00  0.1 ];
    elseif strcmpi(flag,'600')
        data = [ -1.00  0.6
                 -0.83  2.6
                 -0.67  5.8
                 -0.50  7.2
                 -0.33  7.0
                 -0.17  2.6
                  0.00  0.4 ];
    else
        error('%s: flag has to be "400" or "600" for fig. 12!', ...
            upper(mfilename));
    end
% --- FIG 13 ---
elseif fignum==13
    if nargin>2
        error('%s: fig. 13 has no specific data flag.',upper(mfilename));
    else
        % data format: x-axis, -3dB, 3dB, 9dB, 15dB, 25dB
        data = [ -1.0 -10  18  29  37  38
                 -0.9 -10   9  25  31  38
                 -0.7 -10  -2  11  27  33
                 -0.5 -10  -9   4  15  30
                 -0.3 -10  -7   3  13  29
                 -0.1  -7  -1   6  14  29
                  0.1   1   8  13  19  30
                  0.3  10  14  19  26  32
                  0.5  14  20  26  31  34
                  0.7  15  24  30  35  36
                  0.9  11  25  31  34  39
                  1.0   2  19  30  37  38 ];
    end
% --- FIG16 ---
elseif fignum==16
    if nargin>2
        error('%s: fig. 16 has no specific data flag.',upper(mfilename));
    else
        data = [ -1.000  0.0
                 -0.875  3.3
                 -0.750  4.9
                 -0.625  6.6
                 -0.500  7.2
                 -0.375  7.5
                 -0.250  9.0
                 -0.125  9.5
                  0.000  11.1 ];
    end
% --- FIG17 ---
elseif fignum==17
    if nargin>2
        error('%s: fig. 17 has no specific data flag.',upper(mfilename));
    else
        data = [ -9 -0.18 -0.09  0.00  0.09
                 -6 -0.14 -0.05  0.05  0.14
                 -3 -0.09  0.00  0.10  0.18
                  0  0.00  0.09  0.17  0.27
                  3  0.08  0.17  0.23  0.35
                  6  0.10  0.23  0.28  0.37
                  9  0.17  0.25  0.34  0.42 ];
    end
else
    error('%s: no data are avaiable for fig. %i!',upper(mfilename),fignum);
end

