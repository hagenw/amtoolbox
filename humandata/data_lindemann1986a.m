function data = data_lindemann1986a(varargin)
%DATA_LINDEMANN1986a Returns data points from the Lindemann paper
%   Usage: data = data_lindemann1986a_data(flag)
%
%   Output parameters:
%       data    - the data points from the given figure as a matrix with x and y
%                 data or (in cases where no clear x values were given such as
%                 T/2) as a vector containing only the y data
%
%   DATA_LINDEMANN1986A(flag) returns data points from the Lindemann
%   1986a paper. The flag may be one of:
%
%-    fig11_yost - Return data from Fig. 11. with condition 'yost'. The data is
%                  ILD | lateral displacement (10 is max
%                  displacement). If no flag is given, this is the
%                  default.
%
%-    fig11_sayers - Return data from Fig. 11. with condition 'sayers'.
%
%-    fig12_400  - Return data from Fig. 12. where XXX is 400.
%
%-    fig12_600  - Return data from Fig. 12. where XXX is 600.
%
%-    fig13      - Return data from Fig. 13. The output data format is
%                  x-axis, -3dB, 3dB, 9dB, 15dB, 25dB
%
%-    fig16      - Return data from Fig. 16.
%
%-    fig17      - Return data from Fig. 17. XXX what is the data format?
%
%R  lindemann1986a

%   AUTHOR: Hagen Wierstorf


definput.flags.type={...
    'fig11_yost','fig11_sayers',...
    'fig12_400','fig12_600',...
    'fig13',...
    'fig16',...
    'fig17',...
    };
    
  
[flags,keyvals]  = ltfatarghelper({},definput,varargin);


%% ------ Data points from the paper ------------------------------------
% The following data points are guessed after the plots in the paper from
% Lindemann 1986a.
%
% Data for the given figure

if flags.do_fig11_yost
  data = [ 0  0.0
           3  2.1
           6  4.4
           9  6.4
           12  8.6
           15  9.2
           18 10.0 ];
end;

if flags.do_fig11_sayers
  data = [ 6  3.9
           9  4.9
           12  8.9 ];
end;

if flags.do_fig12_400
  data = [ -1.00  0.4
           -0.83  3.3
           -0.67  6.4
           -0.50  8.3
           -0.33  6.7
           -0.17  3.3
           0.00  0.1 ];
end;

if flags.do_fig12_600
  data = [ -1.00  0.6
           -0.83  2.6
           -0.67  5.8
           -0.50  7.2
           -0.33  7.0
           -0.17  2.6
           0.00  0.4 ];
end

if flags.do_fig13
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

if flags.do_fig16
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

if flags.do_fig17
  data = [ -9 -0.18 -0.09  0.00  0.09
           -6 -0.14 -0.05  0.05  0.14
           -3 -0.09  0.00  0.10  0.18
           0  0.00  0.09  0.17  0.27
           3  0.08  0.17  0.23  0.35
           6  0.10  0.23  0.28  0.37
           9  0.17  0.25  0.34  0.42 ];
end

