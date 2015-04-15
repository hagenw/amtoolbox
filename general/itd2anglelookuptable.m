function lookup = itd2anglelookuptable(hrtf,varargin)
%ITD2ANGLELOOKUPTABLE generates an ITD-azimuth lookup table for the given HRTF set
%   Usage: lookup = itd2anglelookuptable(hrtf,fs,model);
%          lookup = itd2anglelookuptable(hrtf,fs);
%          lookup = itd2anglelookuptable(hrtf);
%
%   Input parameters:
%       hrtf   : HRTF data set (as SOFA file or struct)
%       fs     : sampling rate, (default: 44100) / Hz
%       model  : binaural model to use:
%                   'dietz2011' uses the Dietz binaural model (default)
%                   'lindemann1986' uses the Lindemann binaural model
%
%   Output parameters:
%       lookup : struct containing the polinomial fitting data for the
%                ITD -> azimuth transformation, p,MU,S, see help polyfit
%
%   `itd2anglelookuptable(hrtf)` creates a lookup table from the given HRTF data
%   set. This lookup table can be used by the dietz2011 or lindemann1986 binaural
%   models to predict the perceived direction of arrival of an auditory event.
%   The azimuth angle is stored in degree in the lookup table.
%
%   For the handling of the HRTF SOFA file format see
%   http://www.sofaconventions.org/
%
%   See also: dietz2011, lindemann1986, wierstorf2013
%
%   References: dietz2011auditory wierstorf2013 wierstorf2011hrtf

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ===================================
nargmin = 1;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));

definput.flags.model = {'dietz2011','lindemann1986'};
definput.keyvals.fs = 44100;
[flags,kv]=ltfatarghelper({'fs'},definput,varargin);


%% ===== Configuration ==================================================
% Samplingrate
fs = kv.fs;
% Time of noise used for the calculation (samples)
nsamples = fs;
% Noise type to use
noise_type = 'white';
% SFS Toolbox settings
conf.ir.useinterpolation = true;
conf.ir.useoriglength = true;
conf.fs = fs;
conf.c = 343;
conf.usefracdelay = false;


%% ===== Calculation ====================================================
% Generate noise signal
sig_noise = noise(nsamples,1,noise_type);
% Get distance of HRTF set (x0(:,3))
header = sofa_get_header(hrtf);
x0 = SOFAcalculateAPV(header); % x0 = [phi theta r]
% Use azimuth angles in the range from -90deg..90deg
azimuth = rad(-90:90);
nangles = length(azimuth);
% create an empty mod_itd, because the lindemann model didn't use it
mod_itd = [];

if flags.do_dietz2011

    itd = zeros(nangles,12);
    mod_itd = zeros(nangles,23);
    ild = zeros(nangles,23);
    for ii=1:nangles
        % Generate noise coming from the given direction
        % NOTE: this assumes that the ListenerView vector of your HRTF data set
        % is [1 0 0]!
        ir = get_ir(hrtf,[0 0 0],[0 0],[azimuth(ii) 0 x0(ii,3)], ...
                    'spherical',conf);
        sig = auralize_ir(ir,sig_noise,1,conf);
        % Calculate binaural parameters
        [fine, cfreqs, ild_tmp, env] = dietz2011(sig,fs);
        % Unwrap ITD
        itd_tmp = dietz2011unwrapitd(fine.itd,ild_tmp(:,1:12),fine.f_inst,2.5);
        env_itd_tmp = dietz2011unwrapitd(env.itd,ild_tmp(:,13:23),env.f_inst,2.5);
        % Calculate the mean about time of the binaural parameters and store
        % them
        itd(ii,1:12) = median(itd_tmp,1);
        itd(ii,13:23) = median(env_itd_tmp,1);
        ild(ii,:) = median(ild_tmp,1);
    end

elseif flags.do_lindemann1986

    itd = zeros(nangles,36);
    ild = zeros(nangles,36);
    for ii = 1:nangles
        % Generate noise coming from the given direction
        % NOTE: this assumes that the ListenerView vector of your HRTF data set
        % is [1 0 0]!
        ir = get_ir(hrtf,[0 0 0],[0 0],[azimuth(ii) 0 x0(ii,3)], ...
                    'spherical',conf);
        sig = auralize_ir(ir,sig_noise,1,conf);
        % Ten fold upsampling to have a smoother output
        %sig = resample(sig,10*fs,fs);
        % Calculate binaural parameters
        c_s = 0.3; % stationary inhibition
        w_f = 0; % monaural sensitivity
        M_f = 6; % decrease of monaural sensitivity
        T_int = inf; % integration time
        N_1 = 17640; % sample at which first cross-correlation is calculated
        [cc_tmp,dummy,ild(ii,:),cfreqs] = lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1);
        clear dummy;
        cc_tmp = squeeze(cc_tmp);
        % Calculate tau (delay line time) axes
        tau = linspace(-1,1,size(cc_tmp,1));
        % find max in cc
        for jj = 1:size(cc_tmp,2)
            [v,idx] = max(cc_tmp(:,jj));
            itd(ii,jj) = tau(idx)/1000;
        end
    end

end

% Fit the lookup data
for n = 1:size(itd,2)
    [p(:,n),S{n},MU(:,n)] = polyfit(itd(:,n),deg(azimuth)',12);
    [p_ild(:,n),S_ild{n},MU_ild(:,n)] = ...
        polyfit(ild(:,n),deg(azimuth)',12);
end
% Create lookup struct
lookup.p = p;
lookup.MU = MU;
lookup.S = S;
lookup.p_ild = p_ild;
lookup.MU_ild = MU_ild;
lookup.S_ild = S_ild;
