function spatialized = spatialize(input_file,MAP2)
%__________________________________________________________________________
%
% spatialize(catalog,MAP2)
%
% ...
%
% Options
%   input_file      - The HRTF catalog that should be used
%   MAP2            - This a vector in the form [e;a], where a is the
%                     azimuth and e the elevation angle
%
% Outputs:
%   spatialized     - Struct containing the following values:
%                       signal
%                       fs
%                       angle_map
%                       number_of_channels
%                       number_of_directions
%                       number_of_lines (samples?)
%                       
% This functions looks for the used angles (stored in MAP2) and stores the
% corresponding HRTF signals in output_file.
%
% see also: filter_bank, hrtf_analysis
%__________________________________________________________________________
%bi.mo                                                                v.0.3

%
% Jonas Braasch, Wolfgang Heß
% Institut fuer Kommunikationsakustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% e-mail: braasch@ika.ruhr-uni-bochum.de
% 21.07.00
%
% Hagen Wierstorf
% T-Labs, Berlin
% hagen.wierstorf@telekom.de
% 2009/05/28
%

%% Make global variables/constants avaiable
%global def

% disp('-> Loading HRTF catalog');
load(input_file);

% The following two lines make no sense!!!
Fs_cat=Fs;
if Fs~=Fs_cat
   error('sampling frequency of catalog does not match sampling frequency of sound file');
end % of if 


%% Check arguments
if nargin<2
    signal_spatialized = HRIR;
else
    % Look for the right angle entries in the HRTF catalog
    for n = 1:length(MAP2(1,:))
        i = find(MAP(1,:)==MAP2(1,n) & MAP(2,:)==MAP2(2,n));
        signal_spatialized(:,n*2-1:n*2) = HRIR(:,i*2-1:i*2);
    end % of for 
end % of if 

[N_LINES,N_CHANNELS]=size(signal_spatialized);
N_DIR=N_CHANNELS/2;
MAP=MAP2;


%% Returning results
spatialized.signal = signal_spatialized;
spatialized.fs = Fs;
spatialized.angle_map = MAP;
spatialized.number_of_channels = N_CHANNELS;
spatialized.number_of_directions = N_DIR;
spatialized.number_of_lines = N_LINES;

%save([output_file '.sct'],'signal_spatialized','Fs','MAP','N_CHANNELS','N_DIR','N_LINES','-mat');


%% EOF