function [phi,phi_std,itd,ild,cfreqs] = estimate_azimuth(insig,lookup,varargin)
%ESTIMATE_AZIMUTH Estimate the perceived azimuth using a binaural model
%   Usage: [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup,model,do_spectral_weighting,fs)
%          [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup,model,do_spectral_weighting)
%          [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup,model)
%          [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup)
%
%   Input parameters:
%       sig                   : binaural singal
%       lookup                : lookup table to map ITDs to angles (struct)
%       model                 : model to use:
%                                   'dietz2011' (default)
%                                   'lindemann1986'
%       do_spectral_weighting : apply spectral weighting of ITD values after
%                               Raatgever (1980) (default: false)
%       fs                    : sampling rate (default: 44100) (Hz)
%
%   Output parameters:
%       phi     : estimated azimuth (rad)
%       itd     : calculated ITD (s)
%       ild     : calculated ILD (dB)
%       cfreqs  : center frequencies of used auditory filters (Hz)
%
%   `estimate_azimuth(sig,lookup,model,do_spectral_weighting,fs)` uses a
%   binaural model to estimate the perceived direction for a given binaural
%   signal.  Therefore, it needs the struct lookup, which maps ITD values to
%   the corresponding angles. This can be created with the
%   |itd2anglelookuptable| function.  If do_spectral_weighting is set to true,
%   a spectral weighting of the single ITD values after Raatgever is applied. He
%   has done some measurements to see what is the spectral domincance region for
%   lateralization by the ITD and found a region around 600 Hz. Stern et al.
%   have fitted his data with a formula used in this function.
%
%   See also: itd2anglelookuptable
%
%   References: raatgever1980 stern1988 dietz2011auditory lindemann1986a wierstorf2013

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end

if ~isstruct(lookup)
    error('%s: lookup has to be a struct!',upper(mfilename));
end

definput.keyvals.fs = 44100;
definput.flags.binaural_model = {'dietz2011','lindemann1986'};
definput.flags.spectral_weighting = {'no_spectral_weighting','rms_weighting','raatgever_weighting'};
definput.flags.outlier = {'remove_outlier','include_outlier'};

[flags,kv]  = ltfatarghelper({},definput,varargin);


%% ===== Computation ====================================================
%
% === Calculate azimuth values for every frequency channel ===
if flags.do_dietz2011
    ic_threshold=0.98;
    % Run the Dietz model on signal
    [fine,cfreqs,ild] = dietz2011(insig,kv.fs,'nolowpass','fhigh',1400);
    % Unwrap ITDs and get the azimuth values
    itd = dietz2011unwrapitd(fine.itd,ild,fine.f_inst,2.5);
    phi = itd2angle(itd,lookup);
    % Calculate the median over time for every frequency channel of the azimuth
    for n = 1:size(phi,2)
        idx = fine.ic(:,n)>ic_threshold&[diff(fine.ic(:,n))>0; 0]; % compare eq. 9 in Dietz (2011)
        angle = phi(idx,n);
        idx = ~isnan(angle);
        if size(angle(idx),1)==0
            azimuth(n) = NaN;
            azimuth_std(n) = NaN;
        else
            azimuth(n) = median(angle(idx));
            azimuth_std(n) = std(angle(idx));
        end
    end
    % Calculate ITD and ILD values
    itd = median(itd,1);
    % weights for rms-weighting
    rms_weights = fine.rms;
elseif flags.do_lindemann1986
    % run Lindemann model on signal
    c_s = 0.3; % stationary inhibition
    w_f = 0; % monaural sensitivity
    M_f = 6; % decrease of monaural sensitivity
    T_int = inf; % integration time
    N_1 = 1764; % sample at which first cross-correlation is calculated
    [cc_tmp,~,ild,cfreqs] = lindemann1986(insig,kv.fs,c_s,w_f,M_f,T_int,N_1);
    cc_tmp = squeeze(cc_tmp);
    % Calculate tau (delay line time) axes
    tau = linspace(-1,1,size(cc_tmp,1));
    % find max in cc
    itd = zeros(1,size(cc_tmp,2));
    for jj = 1:size(cc_tmp,2)
        [~,idx] = max(cc_tmp(:,jj));
        itd(jj) = tau(idx)/1000;
    end
    azimuth = itd2angle(itd,lookup);
    % TODO: rms weights
end

% === Weights for cross-frequency integration ===
if flags.do_no_spectral_weighting
    w = ones(1,12);
end
if flags.do_rms_weighting
    w = rms_weights;
end
if flags.do_raatgever_weighting
    % Calculate a spectral weighting after Stern1988, after the data of
    % Raatgever1980
    b1 = -9.383e-2;
    b2 =  1.126e-4;
    b3 = -3.992e-8;
    w = 10.^( -(b1*f+b2*(f).^2+b3*(f).^3)/10 );
end

if flags.do_remove_outlier
    % Remove outliers
    [azimuth,azimuth_std,itd,cfreqs,w] = remove_outlier(azimuth,azimuth_std,itd,cfreqs,w);
    %[azimuth_env,~,~] = remove_outlier(azimuth_itd(:,13:23),itd,cfreqs,w);
end

% Calculate mean about frequency channels
if length(azimuth)==0
    phi = NaN;
    phi_std = NaN;
else
    % TODO: add weighted median
    phi = median(azimuth);
    %phi = sum(azimuth.*w)/sum(w);
    %phi_std = std(azimuth);
    phi_std = median(azimuth_std);
end

end % of main function

%% ===== Subfunctions ====================================================
function [azimuth,azimuth_std,itd,cfreqs,w] = remove_outlier(azimuth,azimuth_std,itd,cfreqs,w)
    % remove unvalid ITDs
    %azimuth = azimuth(abs(itd)<0.001);
    %azimuth_std = azimuth_std(abs(itd)<0.001);
    %cfreqs = cfreqs(abs(itd)<0.001);
    %w = w(abs(itd)<0.001);
    %itd = itd(abs(itd)<0.001);
    % remove NaN
    w = w(~isnan(azimuth));
    itd = itd(~isnan(azimuth));
    cfreqs = cfreqs(~isnan(azimuth));
    azimuth = azimuth(~isnan(azimuth));
    azimuth_std = azimuth_std(~isnan(azimuth_std));
    % remove outliers more than 30deg away from median
    if length(azimuth)>0
        itd = itd(abs(azimuth-median(azimuth))<30);
        cfreqs = cfreqs(abs(azimuth-median(azimuth))<30);
        w = w(abs(azimuth-median(azimuth))<30);
        azimuth_std = azimuth_std(abs(azimuth-median(azimuth))<30);
        azimuth = azimuth(abs(azimuth-median(azimuth))<30);
    end
end
