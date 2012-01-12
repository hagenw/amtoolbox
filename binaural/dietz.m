function [hairc_fine_ipd_smooth hairc_mod_ipd_smooth hairc_fine_itd_smooth hairc_mod_itd_smooth hairc_ild hairc_fine_ic hairc_mod_ic hairc_fine_f_inst hairc_mod_f_inst cfreqs] = dietz(signal,fs)
% DIETZ Calculates interaural phase, time and level differences of fine-
%       structure and envelope of the signal, as well as the interaural
%       coherence, which can be used as a weighting function.
%
%   Usage: [fine_ipd,mod_ipd,fine_itd,mod_itd,ild,fine_ic,mod_ic,fine_f_inst,mod_f_inst] =
%          dietz(insig,fs);
%
%   Input parameters:
%       insig       - binaural signal for which values should be calculated
%       fs          - sampling rate (Hz)
%
%   Output parameters:
%
%
%       fine_ipd    - IPD in fine-structure channels
%       mod_ipd     - IPD in envelope channels
%       fine_itd    - ITD in fine-structure channels
%       mod_itd     - ITD in modulation channels
%       ild         - Interaural level difference
%       fine_ic     - Interaural coherence in fine-structure channels
%       mod_ic      - Interaural coherence in modulation channels
%       fine_f_inst - Instantaneous frequency in fine-structure channels
%       mod_f_inst  - Instantaneous frequency in modulation channels
%
%   The steps of the binaural model to calculate the result are the
%   following (see also Dietz et al., 2011):
%
%     1) Middle ear filtering (500-2000 Hz 1st order bandpass)
%
%     2) Auditory bandpass filtering on the basilar membrane using a
%        4th-order all-pole gammatone filterbank, employing 23 filter
%        bands between 200 and 5000 Hz, with a 1 ERB spacing. The filter
%        width was set to correspond to 1 ERB.
%
%     3) Cochlear compression was simulated by power-law compression with
%        an exponent of 0.4.
%
%     4) The transduction process in the inner hair cells was modelled
%        using half-wave rectification followed by filtering with a 770-Hz
%        5th order lowpass.
%
%   The interaural temporal disparities are then extracted using a
%   second-order complex gammatone bandpass (see paper for details).
%
%   Demos: demo_dietz
%

% AUTHOR: Mathias Dietz, Martin Klein-Hennig (implementation for AMT)

% Gammatone filterbank parameters
lower_cutoff_frequency_hz = 200;
upper_cutoff_frequency_hz = 5000;
base_frequency_hz = 1000;
filters_per_ERB = 1;

% Preprocessing parameters
middle_ear_thr = [500 2000]; % Bandpass freqencies for middle ear transfer
middle_ear_order = 2;        % Only even numbers possible
haircell_lp_freq = 770;      % Cutoff frequency for haircell lowpass
haircell_lp_order = 5;       % Order of haircell lowpass
compression_power = 0.4;
alpha = 0;                   % Internal noise strength 65dB = 0.0354
int_noise_case = 'randn';    % randn: add random noise with rms = alpha
% mini: set all values < alpha to alpha

% Parameters for filtering the haircell output
filter_order = 2;            % Used for both mod and fine
filter_attenuation_db = 10;  % Used for both mod and fine

% Only for finestructure plugin
fine_filter_finesse = 3;

% Only for envelope plugin
mod_center_frequency_hz = 216;
mod_filter_finesse = 8;

% For ild- or level-plugin
level_filter_cutoff_hz = 30;
level_filter_order = 2;

% Parameters for localization plugin
coh_param.max_abs_itd   = 1e-3;
coh_param.tau_cycles  = 2.5;    % in cycles cycles
coh_param.tau_s         = 10e-3; % in s for faller

% parameters for data display and analysis
signal_level_dB_SPL = 70; % sound pressure level of left channel

% debugging
displ                = 0;   % display current process (0 = do not)

%% Model processing starts here

dt = 1./fs;
len_in_s = length(signal) ./ fs;
t = linspace(0,len_in_s-dt,length(signal))';

% middle ear band pass filtering
if displ disp('band pass filtering of input according to middle ear transfer charact. -> signal_me'); end
signal_me = middle_ear(signal, middle_ear_thr, middle_ear_order, fs);

% analyzing signals (split into frequency channels)
if displ disp('splitting signal into frequency channels -> signal_filtered'); end

% create filterbank
analyzer = Gfb_Analyzer_new(fs, lower_cutoff_frequency_hz, ...
    base_frequency_hz, upper_cutoff_frequency_hz,...
    filters_per_ERB);
channels = length(analyzer.center_frequencies_hz);
cfreqs=analyzer.center_frequencies_hz;
analyzer_sh = analyzer;

% apply filterbank
[signal_filtered, analyzer] = Gfb_Analyzer_process(analyzer, signal_me(:,1));
[signal_sh_filtered, analyzer_sh] = Gfb_Analyzer_process(analyzer_sh, signal_me(:,2));

% get number of channels
channels = length(cfreqs);

% determine lowpass parameter
tau = coh_param.tau_cycles./cfreqs;

% rectification, comression, and lowpass filtering of filtered signals (haircell processing)
if displ disp('haircell processing of frequency bands -> hairc'); end

% haircell processing
hairc = haircell(signal_filtered, haircell_lp_freq,haircell_lp_order,compression_power, fs)';
hairc_sh = haircell(signal_sh_filtered,haircell_lp_freq,haircell_lp_order,compression_power,fs)';
% also haircell processing, but no 770 Hz filter for fine-structure channel
hairc_nolp = haircell(signal_filtered, '','',compression_power, fs)';
hairc_nolp_sh = haircell(signal_sh_filtered, '','',compression_power, fs)';

% adding internal noise
switch int_noise_case
    case 'randn'
        if displ disp('adding internal random noise -> hairc'); end
        hairc = hairc + (randn(length(hairc),channels)*alpha);
        hairc_sh = hairc_sh + (randn(length(hairc_sh),length(cfreqs))*alpha);
        hairc_nolp = hairc_nolp + (randn(length(hairc),channels)*alpha);
        hairc_nolp_sh = hairc_nolp_sh + (randn(length(hairc_sh),length(cfreqs))*alpha);
    case 'mini'
        if displ disp('adding internal noise via minimum -> hairc'); end
        hairc = max(hairc,alpha);
        hairc_sh = max(hairc_sh,alpha);
        hairc_nolp = max(hairc,alpha);
        hairc_nolp_sh = max(hairc_sh,alpha);
    otherwise
        disp('Unknown noise case. WARNING: No noise added!')
end

% processing the hairc output with a modulation frequency filter
%cmin = min(find(cfreqs>2*mod_center_frequency_hz)); % lowest freq. band for envelope detection
if displ disp('enveloping haircell output -> hairc_mod'); end
mod_filter_bandwidth_hz = mod_center_frequency_hz/mod_filter_finesse;
[hairc_mod, hairc_sh_mod] =... %Gfb_envelope_filter(hairc(:,cmin:end), hairc_sh(:,cmin:end), fs,...
    Gfb_envelope_filter(hairc, hairc_sh, fs,...
    mod_center_frequency_hz, mod_filter_bandwidth_hz, ...
    filter_attenuation_db, filter_order);

% calculation of interaural functions from haircell modulation
if displ disp('calculating interaural functions from haircell modulation'); end
[hairc_mod_itf, hairc_mod_ipd, hairc_mod_ipd_smooth,...
    hairc_mod_ild, hairc_mod_itd, hairc_mod_itd_C, hairc_mod_itd_smooth,...
    hairc_mod_itd_C_smooth, hairc_mod_f_inst,...
    hairc_mod_level, hairc_mod_ic, hairc_mod_rms] =...
    interaural_functions(hairc_mod, hairc_sh_mod, tau,...
    mod_center_frequency_hz+0*cfreqs,...
    signal_level_dB_SPL, compression_power, coh_param.tau_cycles, fs);

% processing the hairc output with a fine structure filter
if displ disp('deriving fine structure of haircell output -> hairc_fine'); end
[hairc_fine, hairc_sh_fine] =...
    Gfb_envelope_filter(hairc_nolp, hairc_nolp_sh, fs, cfreqs,...
    cfreqs/fine_filter_finesse, filter_attenuation_db, filter_order);

% calculation of interaural functions from haircell fine structure
if displ disp('calculating interaural functions from haircell fine structure'); end
[hairc_fine_itf, hairc_fine_ipd, hairc_fine_ipd_smooth,...
    hairc_fine_ild, hairc_fine_itd, hairc_fine_itd_C, hairc_fine_itd_smooth,...
    hairc_fine_itd_C_smooth, hairc_fine_f_inst,...
    hairc_fine_level, hairc_fine_ic, hairc_fine_rms] =...
    interaural_functions(hairc_fine, hairc_sh_fine, tau, cfreqs,...
    signal_level_dB_SPL, compression_power, coh_param.tau_cycles, fs);

% determine ILD of the hairc output
if displ disp('determining ild of the haircell output -> hairc_ild'); end
[hairc_ild] =...
    ild_filter(hairc,hairc_sh,level_filter_cutoff_hz,...
    level_filter_order,compression_power,fs);

% remove finestructure information > 1400 Hz
hairc_fine_f_inst(:,cfreqs>1400)=[];
hairc_fine_ic(:,cfreqs>1400)=[];
hairc_fine_ipd_smooth(:,cfreqs>1400)=[];
hairc_fine_itd_smooth(:,cfreqs>1400)=[];


if displ disp('finished');end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% middle_ear %%%%%%%%%%%%%%%%%%%%%%%%
function output = middle_ear(input,thresholds,order,fs)

% bandpass filtering
[b,a] = butter(order,thresholds(1)/(fs/2),'low');
low_filtered= filter(b,a,input);

[b,a] = butter(order,thresholds(2)/(fs/2),'high');
output= filter(b,a,low_filtered);

end

%% Gfb_envelope_filter %%%%%%%%%%%%%%%%%
function [envelopes_filtered, envelopes_sh_filtered] = Gfb_envelope_filter(s1, s2, sampling_rate_hz, center_frequency_hz,...
    bandwidth_hz, attenuation_db, gamma_filter_order)
% [envelopes_filtered, envelopes_sh_filtered] =...
%   Gfb_envelope_filter(s1, s2, sampling_rate_hz, center_frequency_hz, bandwidth_hz, attenuation_db, gamma_filter_order);
%
% Filters each row of s1 and s2 with the gammatone filter defined by the input parameters.
% Takes both vectors and matrices.
%
% Input
%   s1, s2 - signals to be filtered
%   sampling_rate_hz - sampling frequency / Hz
%   center_frequency_hz - centre frequency of the gammatone filter / Hz
%   bandwidth_hz - bandwidth of the gammatone filter at the level attenuation_db / Hz
%   attenuation_db - attenuation in dB at which the filter has bandwidth_hz
%   gamma_filter_order - order of the filter

s1 = s1';
s2 = s2';

[M, N] = size(s1);
if length(center_frequency_hz) == 1
    center_frequency_hz = center_frequency_hz * ones(1,M);
end
if isempty(bandwidth_hz) % default: width = 1 ERB
    recip_width1erb = diff(Gfb_hz2erbscale(1:N/2));
    bandwidth_hz = round(1./recip_width1erb(round(center_frequency_hz)));
elseif length(bandwidth_hz) == 1
    bandwidth_hz = bandwidth_hz * ones(1,M);
end

for i = 1:M
    filter = Gfb_Filter_new(sampling_rate_hz, center_frequency_hz(i),...
        bandwidth_hz(i), attenuation_db, gamma_filter_order);
    [envelopes_filtered(:,i),    filter_obj] = Gfb_Filter_process(filter, s1(i,:));
    [envelopes_sh_filtered(:,i), filter_obj] = Gfb_Filter_process(filter, s2(i,:));
end
end

%% interaural_functions
function [itf, ipd, ipd_lp, ild, itd, itd_C, itd_lp, itd_C_lp,...
    finst, level, ic, weight] =...
    interaural_functions(s1, s2, tau, cfreqs, signal_level_dB_SPL,...
    compr, coh_cycles, fs)
% Calculates interaural parameters of two complex-valued
% Hilbert-transformed signals, supplied by the GFB Analyzer()
%
% Input
%   s1, s2 - input signals
%   tau - lowpass filter parameter, such that a = exp(-1./(fs*tau)), lowpass = filter([1-a], [1, -a], x)
%   Note: If tau is a scalar, lowpass is done on every channel with the value of tau.
%         If the filter parameter tau is a vector, it has to have as many elements
%         as the number of channels of s1 and s2; each value tau(i) is applied to
%         the signals in s1(i,:) and s2(i,:).
%   cfreqs - central frequencies
%   fs - sampling frequencies
% Output
%   itf - transfer function
%   itf_equal - transfer function without amplitude
%   ipd - phase difference in rad
%   ipd_lp - based on lowpass-filtered itf, phase difference in rad
%   ild - level difference in dB
%   itd, itd_C, itd_lp, itd_C_lp - time difference based on instantaneous
%                and central frequencies, with and without low-passed itf
%   finst_1 - instantaneous frequencies in the channels of the filtered s1
%   finst_2 - instantaneous frequencies in the channels of the filtered s2

a = exp( -1./(fs*tau) );

itf = s2 .* conj(s1);

ipd_lp = zeros(size(itf));
finst_1 = zeros(size(itf));
finst_2 = zeros(size(itf));

% interaural phase difference
ipd = angle(itf);

% interaural coherence
ic = new_ic(itf, coh_cycles./cfreqs, fs);

% interaural level difference for the case of tau - scalar

s1_low = lowpass(abs(s1),a);
s2_low = lowpass(abs(s2),a);                % take envelope at higher frequencies
s2_low( find(abs(s2_low) < eps) ) = eps;    % avoid division by zero
s1_low( find(abs(s1_low) == 0 ) ) = eps;    % avoid log(0)
ild = 20*log10(s1_low./s2_low)./compr;
if length(a) == 1
    a(1:length(cfreqs)) = a;
    tau(1:length(cfreqs)) = tau;
end
for k = 1:length(cfreqs);
    ipd_lp(:,k) = angle(lowpass(itf(:,k),a(k)))';
end

% interaural time difference, based on central and instantaneous frequencies
for k = 1:length(cfreqs)
    finst_1(:,k) = calc_finst(s1(:,k),fs,tau(k),0);
    finst_2(:,k) = calc_finst(s2(:,k),fs,tau(k),0);
    itd_C(:,k) = 1/(2*pi)*ipd(:,k)/cfreqs(k);
    itd_C_lp(:,k) = 1/(2*pi)*ipd_lp(:,k)/cfreqs(k);
end
finst = max(eps,0.5*(finst_1 + finst_2));

% to avoid division by zero

itd = 1/(2*pi)*ipd./finst;    % based on instantaneous frequencies
itd_lp = 1/(2*pi)*ipd_lp./finst;    % based on instantaneous frequencies

level = signal_level_dB_SPL*compr + 20*log10(abs(s1)+abs(s2));



% weighting of channels for cumulative ixd determination
weight = signal_level_dB_SPL*compr + 20*log10(sqrt(2)*min(rms(abs(s1)),rms(abs(s2))));
% sqrt(2) is due to half-wave rectification (included 28th Sep 07)
weight = max(weight,0); % avoid negative weights
end

%% lowpass
function y = lowpass(x, a)
% y = lowpass(x, a)
% This is a simple low-pass filter y(n) = (1-a)*x(n) - a*y(n-1)
% Meaning of parameter a:
%   a - damping coefficient (0 - no filtering, 1 - flat output)
%   tau = 1/(2*pi*f_c)      where f_c is the cutoff frequency of the filter
%   a = exp(-1/(f_s*tau))   where fs - sampling frequency
%
% Input
%   x - input signal
%   1] The signal x may be either row or column vector, the output y is ALWAYS a column vector
%   2] If x is a matrix, time is considered along the LONGER dimension of x
%      Then, the output y is always a column vector or a matrix with time along the rows
%
% Output
%  y - filtered signal
%
% Example see ...\examples\example_lowpass_tester.m

[rows, columns] = size(x);
if rows < columns
    x = x.';
    y = zeros(columns,rows);
else
    y = zeros(rows, columns);
end
[rows, columns] = size(y);

y = filter([1-a], [1, -a], x);
end

%% finst - instantaneous frequency
function finst = calc_finst(sig,fs,tau,norm)
%
% function finst = calc_finst(sig,fs,tau,norm);
%
% Calculates instantaneous frequency from a complex (analytical) signal
% using first order differences
%
% input parameters:
% sig:  complex (analytical) input signal
% fs:   sampling frequency of sig
% tau:  exponential decay time of temporal averaging filter
% norm: exponent for amplitude weighting for temporal averaging filter
%       process(0: no level weigthing)
%
% output values:
% finst:   vector of estimated inst. frequency values (w temporal
%          averaging
%
% copyright: Universitaet Oldenburg
% author   : vh
% date     : 12/04
%

sig = sig';

alpha = exp(-1/tau/fs);
b = [1-alpha];
a = [1 -alpha];

finst = sig./(abs(sig)+eps);
finst = abs(sig).^norm.*[0 finst(2:end).*conj(finst(1:end-1))];
finst = filter(b,a,finst);
finst = angle(finst')/2/pi*fs;

end

function ic = new_ic(itf,tau_coherence,fs)

%tau_coherence = 15e-3; % good value for ipd_fine

c_coh = exp(-1./(fs.*tau_coherence));
if length(tau_coherence)==1
    ic = abs(filter(1-c_coh,[1 -c_coh],itf))./abs(filter(1-c_coh,[1 -c_coh],abs(itf)));
elseif length(tau_coherence)==size(itf,2)
    ic = zeros(size(itf));
    for n = 1:length(tau_coherence)
        ic(:,n) = abs(filter(1-c_coh(n),[1 -c_coh(n)],itf(:,n)))./ ...
            abs(filter(1-c_coh(n),[1 -c_coh(n)],abs(itf(:,n))));
    end
else
    error('wrong number of tau_coherence values')
end
end

%% ild_filter
function output = ild_filter(hairc,hairc_sh,lp_threshold_freq,lp_order,compression,fs)

% lowpass filtering
[b,a] = butter(lp_order,lp_threshold_freq/(fs/2),'low');
hclp    = filter(b,a,hairc);
hclp_sh = filter(b,a,hairc_sh);

output = 20*log10(max(hclp_sh,1e-4)./max(hclp,1e-4))/compression;
end

%% haircell
function output = haircell(input,cutoff,order,compress_power,fs)

% half-wave rectification
rect = max(real(input),0);

% compression
output = rect.^compress_power;

% lowpass filtering, only if desired
if (cutoff~='' | order ~= '')
    fNorm = cutoff*(1/sqrt((2^(1/order)-1))) / (fs/2);
    [b,a] = butter(1,fNorm,'low');
    for k = 1:order
        output= filter(b,a,output);
    end
end

end
%% rms
function y = rms(x)

y = sqrt(sum(x .^ 2) ./ length(x));
end

%% ------------------------------------------------------------------------
%%-----------------------------------------------------------------------------
%%
%%   Copyright (C) 2002-2011   AG Medizinische Physik,
%%                             Universitaet Oldenburg, Germany
%%                             http://www.physik.uni-oldenburg.de/docs/medi
%%
%%   Permission to use, copy, and distribute this software/file and its
%%   documentation for any purpose without permission by UNIVERSITAET OLDENBURG
%%   is not granted.
%%
%%   Permission to use this software for academic purposes is generally
%%   granted.
%%
%%   Permission to modify the software is granted, but not the right to
%%   distribute the modified code.
%%
%%   This software is provided "as is" without expressed or implied warranty.
%%
%%   Authors: Tobias Peters (tobias@medi.physik.uni-oldenburg.de) 2002
%%            Mathias Dietz (mathias.dietz@uni-oldenburg.de)      2006-2009
%%            Martin Klein-Hennig (martin.klein.hennig@uni-oldenburg.de)
%%            2011
%%
%%
%%-------------------------------------------------------------------------

%OLDFORMAT
