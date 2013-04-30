function lookup = itdazimuthlookuptable(irs,varargin)
%ITDAZIMUTHLOOKUPTABLE generates an ITD-azimuth lookup table for the given HRTF set
%   Usage: lookup = itdazimuthlookuptable(irs,fs,model);
%          lookup = itdazimuthlookuptable(irs,fs);
%          lookup = itdazimuthlookuptable(irs);
%
%   Input parameters:
%       irs    : HRTF data set (TU Berlin irs format)
%       fs     : sampling rate, default: 44100 (Hz)
%       model  : binaural model to use, default: 'dietz'                   
%
%   Output parameters:
%       lookup : struct containing lookup data
%
%   `itdazimuthlookuptable(irs)` creates a lookup table from the given IR data
%   set. This lookup table can be used by the dietz binaural model to predict
%   the perceived direction of arrival of an auditory event.
%   
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ===================================
nargmin = 1;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));

definput.flags.model = {'dietz','lindemann'};
definput.keyvals.fs = 44100;
[flags,kv]=ltfatarghelper({'fs'},definput,varargin);


%% ===== Configuration ==================================================
% time of noise used for the calculation (samples)
nsamples = kv.fs;
% noise type to use
noise_type = 'white';


%% ===== Calculation ====================================================
% generate noise signal
sig_noise = noise(nsamples,1,noise_type);
% get only the -90 to 90 degree part of the irs set
idx = (( irs.apparent_azimuth>-pi/2 & irs.apparent_azimuth<pi/2 & ...
    irs.apparent_elevation==0 ));
irs = slice_irs(irs,idx);
% iterate over azimuth angles
nangles = length(irs.apparent_azimuth);
% create an empty mod_itd, because the lindemann model didn't use it
mod_itd = [];

if flags.do_dietz

    itd = zeros(nangles,12);
    mod_itd = zeros(nangles,23);
    ild = zeros(nangles,23);
    for ii = 1:nangles
        % generate noise coming from the given direction
        ir = get_ir(irs,irs.apparent_azimuth(ii));
        sig = auralize_ir(ir,sig_noise);
        % calculate binaural parameters
        [fine, modulation, cfreqs, ild_tmp] = dietz2011(sig,fs);
        % unwrap ITD
        itd_tmp = unwrap_itd(fine.itd(:,1:12),ild_tmp,fine.f_inst);
        % calculate the mean about time of the binaural parameters and store
        % them
        itd(ii,:) = median(itd_tmp,1);
        mod_itd(ii,:) = median(modulation.itd,1);
        ild(ii,:) = median(ild_tmp,1);
    end

elseif flags.do_lindemann

    itd = zeros(nangles,36);
    ild = zeros(nangles,36);
    for ii = 1:nangles
        % generate noise coming from the given direction
        ir = get_ir(irs,irs.apparent_azimuth(ii));
        sig = auralize_ir(ir,sig_noise);
        % Ten fold upsampling to have a smoother output
        %sig = resample(sig,10*fs,fs);
        % calculate binaural parameters
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

% Create lookup struct
lookup.itd = itd;
lookup.ild = ild;
lookup.mod_itd = mod_itd;
lookup.cfreq = cfreqs;
lookup.azimuth = irs.apparent_azimuth;
