function [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup,model,fs)
%ESTIMATE_AZIMUTH uses a binaural model to estimate the perceived azimuth angle
%
%   Usage: [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup,model,fs)
%          [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup,model)
%          [phi,itd,ild,cfreqs] = estimate_azimuth(sig,lookup)
%
%   Input parameters:
%       sig     - binaural singal
%       lookup  - lookup table to map ITDs to angles (struct)
%       model   - model to use:
%                   'dietz' (default)
%                   'lindemann'
%       fs      - sampling rate (default: 44100) (Hz)
%
%   Output parameters:
%       phi     - estimated azimuth (rad)
%       itd     - calculated ITD (s)
%       ild     - calculated ILD (dB)
%       cfreqs  - center frequencies of used auditory filters (Hz)
%
%   ESTIMATE_AZIMUTH(sig,lookup) uses a binaural model to estimate the perceived
%   direction for a given binaural signal. therefore it needs the struct lookup,
%   which maps ITD values to the corresponding angles.
%
%   see also: lookup

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));

if nargin==2
    model = 'dietz';
    fs = 44100;
elseif nargin==3
    isargstring(model);
    fs = 44100;
else
    isargchar(model);
    isargpositivescalar(fs);
end


%% ===== Computation ====================================================
%
% === Dietz model ===
if strcmpi('dietz',model)

    ic_threshold=0.98;
    cn = 10; % channel number for time plot

    % TODO: the new Matlab Syntax for dummy variables is:
    % [~,idx] = max(a);
    %max This should also work in Octave, but version 4.3 is needed. Check if this
    % could be installed easily under Debian

    % Run the Dietz model on signal
    [fine, modulation, cfreqs, ild_tmp] = dietz2011(sig,fs);

    % Unwrap ITDs and get the azimuth values
    itd = unwrap_itd(fine.itd(:,1:12),ild_tmp,fine.f_inst);
    phi = itd2azimuth(itd,lookup);

    % Calculate the median over time for every frequency channel of the azimuth
    for n = 1:size(phi,2)
        idx = fine.ic(:,n)>ic_threshold&[diff(fine.ic(:,n))>0; 0];
        %phi(:,n) = phi(fine_ic>ic_threshold&[diff(fine_ic)>0; zeros(1,12)],n);
        angle = phi(idx,n);
        idx = ~isnan(angle);
        %azimuth(n) = mean(angle(idx));
        if size(angle(idx),1)==0
            azimuth(n) = NaN;
        else
            azimuth(n) = median(angle(idx));
        end
    end
    % Calculate ITD and ILD values
    ild = median(ild_tmp,1);
    itd = median(itd,1);
    %w = spectral_weight(cfreqs(1:12));
    %for ii=1:size(itd,1)
    %    itd_mean(ii) = sum(w.*itd(ii,:))./sum(w);
    %end
    %phi_mean = itdmean2azimuth(itd_mean,lookup);

elseif strcmpi('lindemann',model)

    % run Lindemann model on signal
    c_s = 0.3; % stationary inhibition
    w_f = 0; % monaural sensitivity
    M_f = 6; % decrease of monaural sensitivity
    T_int = inf; % integration time
    N_1 = 1764; % sample at which first cross-correlation is calculated
    [cc_tmp,dummy,ild,cfreqs] = lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1);
    clear dummy;
    cc_tmp = squeeze(cc_tmp);
    % Calculate tau (delay line time) axes
    tau = linspace(-1,1,size(cc_tmp,1));
    % find max in cc
    itd = zeros(1,size(cc_tmp,2));
    for jj = 1:size(cc_tmp,2)
        [v,idx] = max(cc_tmp(:,jj));
        itd(jj) = tau(idx)/1000;
    end
    azimuth = itd2azimuth(itd,lookup);
end

% NOTE: the last version used the mean instead of median, and no outliers had
% been removed!
%
% Calculate mean about frequency channels
% first remove outliers
[azimuth,cfreqs] = remove_outlier(azimuth,itd,cfreqs);
%w = spectral_weighting(cfreqs);
%phi = sum(azimuth.*w)/sum(w);
%phi = phi_mean;
if length(azimuth)==0
    phi = rad(90);
else
    phi = median(azimuth);
end

end % of main function

%% ===== Subfunctions ====================================================
function [azimuth,cfreqs] = remove_outlier(azimuth,itd,cfreqs)
    % remove unvalid ITDs
    azimuth = azimuth(abs(itd(1:12))<0.001);
    % remove NaN
    azimuth = azimuth(~isnan(azimuth));
    % remove azimuth with other sign than the mean
    %azimuth = azimuth(sign(azimuth)==sign(median(azimuth)));
    % remove outliers with 2*std
    %azimuth = azimuth(abs(azimuth-median(azimuth))<std(azimuth));
    azimuth = azimuth(abs(azimuth-median(azimuth))<rad(30));
    %idx = (( abs(itd(1:12))<0.001 && ...
    %         ~isnan(azimuth) && ...
    %         abs(azimuth-median(azimuth))<rad(30) ));
    %azimuth = azimuth(idx);
    %cfreqs = cfreqs(idx);
end
function w = spectral_weighting(f)
    % Calculate a spectral weighting after Stern1988, after the data of
    % Raatgever1980
    b1 = -9.383e-2;
    b2 =  1.126e-4;
    b3 = -3.992e-8;
    w = 10.^( -(b1*f+b2*(f).^2+b3*(f).^3)/10 );
end
