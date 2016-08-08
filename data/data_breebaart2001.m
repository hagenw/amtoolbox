function data = data_breebaart2001(varargin)
%DATA_BREEBAART2001 Returns data points from the Breebaart et al. (1999) paper
%   Usage: data = data_breebaart(figure,nfc)
%
%   The figure may be one of:
%   'fig3'     (default) Returns the N0Spi values for figure 3
%   'fig6'      Returns the NpiS0 values for figure 6
%
%   The nfc (center frequency of noise) may be one of: 
%   'nfc125','nfc250','nfc500','nfc1000'    -> fig3 and fig6
%   'nfc2000','nfc4000' -> fig3 only         
%
%   Examples:
%   ---------
% 
%   To get data for the fig. 3 Breebaart et al. for the 
%   condition with 125 Hz center frequency use :::
%
%     data_breebaart2001('fig3','nfc125');
%
%   References: 


%   AUTHOR: Martina Kreuzbichler


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type={'fig3','fig6'};
definput.flags.nfc = {'nfc125','nfc250','nfc500','nfc1000','nfc2000','nfc4000'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%% ------ Data points from the paper ------------------------------------
%
% Data for the given figure
if flags.do_fig3
    if flags.do_nfc125
        data= [-20.5 -21 -21.5 -22 -21 -23.5];
    elseif flags.do_nfc250
        data = [-22 -22 -20.5 -21 -21 -24 -30];
    elseif flags.do_nfc500
        data = [-20 -21 -21 -20 -21 -22 -26 -28.5];
    elseif flags.do_nfc1000
        data = [-18 -17.5 -18 -17 -18 -18 -20 -22 -25.5];
    elseif flags.do_nfc2000
        data = [-13.5 -10.5 -11 -14 -14 -15 -15 -17.5 -18 -23];
    elseif flags.do_nfc400
        data = [-9.5 -10 -11.5 -14 -15.5 -17 -16.5 -16 -17.5 -21];
    end
elseif flags.do_fig6
    if flags.do_nfc125
        data= [-12.5 -12 -11.5 -8.5 -6 -8.5];
    elseif flags.do_nfc250
        data = [-17 -16.5 -18.5 -19 -16 -20 -24];
    elseif flags.do_nfc500
        data = [-16.5 -19 -19 -20 -17.5 -18.5 -21.5 -25.5];
    elseif flags.do_nfc1000
        data = [-15 -17.5 -14.5 -16 -16.5 -15.5 -17 -21.5 -23.5];
    end
    
end
    
    
