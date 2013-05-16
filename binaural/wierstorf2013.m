function [localization_error,perceived_direction,desired_direction,x0] = ...
        wierstorf2013(X,phi,xs,src,L,method,number_of_speakers,array);
%WIERSTORF2013 estimate the localization within a WFS or stereo setup
%   Usage: [loc_error,...] = wierstorf2013(X,phi,xs,src,L,'wfs',nls,array);
%          [loc_error,...] = wierstorf2013(X,phi,xs,src,L,'stereo');
%
%   Input parameters:
%       X           : position of the listener [x,y] in m
%       phi         : orientation of the listener in rad (0 is in the direction
%                     of the x-axis
%       xs          : position of the point source in m / direction of the 
%                     plane wave
%       src         : source type
%                       'ps' for a point source
%                       'pw' for a plane wave
%       L           : length/diameter of the loudspeaker array (WFS only)
%       method      : reproduction setup
%                       'wfs' for wave field synthesis
%                       'setreo' for stereophony
%       nls         : number of loudspeakers of the array (WFS only)
%       array       : loudspeaker array type (WFS only);
%                       'linear' for a linear loudspeaker array
%                       'circle' for a circular loudspeaker array
%
%   Output parameters:
%       localization_error  : deviation from the desired direction, defined as
%                             perceived_direction - desired_direction / rad
%       perceived_direction : the direction of arrival the binaural model has
%                             estimated for our given setup / rad
%       desired_direction   : the desired direction of arrival indicated by the
%                             source position xs / rad
%       x0                  : position and directions of the loudspeakers in the
%                             form n x 6, where n is the number of loudspeakers
%
%   `wierstorf2013(X,phi,xs,src,'wfs',L,nls,array)` calculates the localization
%   error for the defined wave field synthesis or stereophony setup. The
%   localization error is defined here as the difference between the perceived
%   direction as predicted by the dietz2011 binaural model and the desired
%   direction given by *xs*. The loudspeaker setup for the desired reproduction
%   method is simulated via HRTFs which are than convolved with white noise
%   which is fed into the binaural model.
%
%   For the simulation of the wave field synthesis or stereophony setup this
%   functions depends on the Sound-Field-Synthesis Toolbox, which is available
%   here: http://github.com/sfstoolbox/sfs. It runs under Matlab and Octave. The
%   revision used to genrate the figures in the corressponding paper is
%   a8914700a4.
%
%   See also: estimate_azimuth, dietz2011
%
%   References: wierstorf2011hrtf wierstorf2013

% AUTHOR: Hagen Wierstorf

%   Copyright (c) 2013   Assessment of IP-based Applications
%                        Technische Universitaet Berlin
%                        Ernst-Reuter-Platz 7, 10587 Berlin, Germany


%% ===== Checking of input parameters and dependencies ===================
nargmin = 6;
nargmax = 8;
narginchk(nargmin,nargmax);

% Checking for the Sound-Field-Synthesis Toolbox
if !which('SFS_start')
    error(['%s: you need to install the Sound-Field-Synthesis Toolbox.\n', ...
        'You can download it at https://github.com/sfstoolbox/sfs.\n', ...
        'The results in the paper are verified up to revision a8914700a4'], ...
        upper(mfilename));
end


%% ===== Configuration ===================================================
fs = 44100; % / Hz

%% ===== Loading of additional data ======================================
% load HRTFs, see:
% https://dev.qu.tu-berlin.de/projects/measurements/wiki/2010-11-kemar-anechoic
[~,path] = download_hrtf('wierstorf2011_3m');
load([path 'wierstorf2011_3m.mat']);
hrtf = irs;
% load lookup table to map ITD values of the model to azimuth angles.
% the lookup table was created using the same HRTF database
load([path(1:end-10) 'modelstages/wierstorf2013itd2anglelookup.mat']);


%% ===== Simulate the binaural ear signals ===============================
% Setting up the configuration
conf = SFS_config;
% Simulate a stereo setup
if strcmpi('stereo',method)
    number_of_speakers = 2;
    array = 'linear';
elseif !strcmpi('wfs',method)
    error('%s: %s is not a valid method.',upper(mfilename),method);
end
conf.array = array;
% calculating the distance between the loudspeakers
if strcmpi('circle',array)
    conf.dx0 = pi*L/number_of_speakers;
elseif strcmpi('linear',array)
    conf.dx0 = L/(number_of_speakers-1);
else
    error('%s: %s is not a valid array.',upper(mfilename),array);
end
% get loudspeaker positions
x0 = secondary_source_positions(L,conf);
% selection of loudspeakers for WFS
if strcmpi('wfs',method) && strcmpi('circle',array)
    x0 = secondary_source_selection(x0,xs,src);
end
% simulate the binaural impulse response
if strcmpi('stereo',method)
    % first loudspeaker
    ir1 = ir_point_source(X,phi,x0(1,1:3),hrtf,conf);
    % second loudspeaker
    ir2 = ir_point_source(X,phi,x0(2,1:3),hrtf,conf);
    % sum of both loudspeakers
    ir = (ir1+ir2)/2;
else % WFS
    ir = ir_wfs_25d(X,phi,xs,src,L,hrtf,conf);
end
% FIXME: this should have no influence on the results
% scale impulse response
%ir = 0.9 * (ir ./ (max(abs(ir(:)))+eps));
% generate a 0.1s noise signal
sig_noise = noise(fs/10,1,'white');
% convolve with impulse response
sig = auralize_ir(ir,sig_noise);


%% ===== Estimate the direction of arrival for the listener ==============
% this is done by calculating ITDs with the dietz2011 binaural model, which are
% then mapped to azimuth values with a lookup table
%
% estimate the perceived direction of arrival
perceived_direction = estimate_azimuth(sig,lookup,'dietz2011',0);
% calculate the desired direction
desired_direction = source_direction(X,phi,xs,src);
% calculate the localization error as the difference of both
localization_error = perceived_direction - desired_direction;

end % of main function


%% ----- Subfunctions ----------------------------------------------------
function direction = source_direction(X,phi,xs,src)
    if strcmp('pw',src)
        [direction,~,~] = cart2sph(xs(1),xs(2),0);
    elseif strcmp('ps',src)
        x = xs-X;
        [tmp,~,~] = cart2sph(x(1),x(2),0);
        direction = tmp-pi;
    end
    direction = direction + phi;
end
