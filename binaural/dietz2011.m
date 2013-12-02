function [fine, env, fc, ild] = dietz2011(insig,fs,varargin)
%DIETZ2011  Dietz 2011 binaural model
%   Usage: [...] = dietz(insig,fs);
%
%   Input parameters:
%       insig       : binaural signal for which values should be calculated
%       fs          : sampling rate (Hz)
%
%   Output parameters:
%       fine        : Information about the fine structure (see below)
%       env         : Information about the envelope (see below)
%
%   `dietz2011(insig,fs)` calculates interaural phase, time and level
%   differences of fine- structure and envelope of the signal, as well as
%   the interaural coherence, which can be used as a weighting function.
%
%   The output structures *fine* and *env* have the following fields:
%
%     s1        : left signal as put in the binaural processor
%     s2        : right signal as put in the binaural processor
%     fc        : center frequencies of the channels (f_carrier or f_mod)
%     itf       : transfer function
%     itf_equal : transfer function without amplitude
%     ipd : phase difference in rad
%     ipd_lp    : based on lowpass-filtered itf, phase difference in rad
%     ild : level difference in dB
%     itd, itd_C, itd_lp, itd_C_lp - time difference based on instantaneous
%                  and central frequencies, with and without low-passed itf
%     f_inst_1 : instantaneous frequencies in the channels of the filtered s1
%     f_inst_2 : instantaneous frequencies in the channels of the filtered s2
%     f_inst   : instantaneous frequencies (average of f_inst1 and 2)
%  
%   The steps of the binaural model to calculate the result are the
%   following (see also Dietz et al., 2011):
%
%   1) Middle ear filtering (500-2000 Hz 1st order bandpass)
%
%   2) Auditory bandpass filtering on the basilar membrane using a
%      4th-order all-pole gammatone filterbank, employing 23 filter
%      bands between 200 and 5000 Hz, with a 1 ERB spacing. The filter
%      width was set to correspond to 1 ERB.
%
%   3) Cochlear compression was simulated by power-law compression with
%      an exponent of 0.4.
%
%   4) The transduction process in the inner hair cells was modelled
%      using half-wave rectification followed by filtering with a 770-Hz
%      5th order lowpass.
%
%   The interaural temporal disparities are then extracted using a
%   second-order complex gammatone bandpass (see paper for details).
%
%   `dietz2011` accepts the following optional parameters:
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 200 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 5000 Hz.
%
%     'basef',basef  Ensure that the frequency basef is a center frequency
%                    in the filterbank. The default value  is 1000.
%
%     'filters_per_ERB', filters_per_erb
%                    Filters per erb. The default value is 1.
%
%     'middle_ear_thr',r
%                    Bandpass freqencies for middle ear transfer. The
%                    default value is `[500 2000]`.
%
%     'middle_ear_order',n 
%                    Order of middle ear filter. Only even numbers are
%                    possible. The default value is 2.
%
%     'compression_power',cpwr
%                    XXX. The default value is 0.4.
%
%     'alpha',alpha  Internal noise strength. Convention FIXME 65dB =
%                    0.0354. The default value is 0.
%
%     'int_randn'    Internal noise XXX. This is the default.
%
%     'int_mini'     Internal noise XXX.
%
%     'filter_order',fo
%                    Filter order for output XXX. Used for both 'mod' and 'fine'.
%                    The default value is 2.
%
%     'filter_attenuation_db',fadb
%                    XXX. Used for both 'mod' and 'fine'. The default value is 10.
% 
%     'fine_filter_finesse',fff
%                    Only for finestructure plugin. The default value is 3.
%
%     'mod_center_frequency_hz',mcf_hz
%                    XXX. Only for envelope plugin. The default value is 135.
%
%     'mod_filter_finesse',mff
%                    XXX. Only for envelope plugin. The default value is 8.
% 
%     'level_filter_cutoff_hz',lfc_hz
%                    XXX. For ild- or level-plugin. The default value is 30.
%
%     'level_filter_order',lforder
%                    XXX. For ild- or level-plugin. The default value is 2.
%
%     'coh_param',coh_param
%                    This is a structure used for the localization
%                    plugin. It has the following fields:
%  
%                      `max_abs_itd`
%                         XXX. The default value is 1e-3.
%
%                      `tau_cycles`
%                         XXX. The default value is 5. 
%
%                      `tau_s`
%                         XXX. The default value is 10e-3.
%
%     'signal_level_dB_SPL',signal_level
%                    Sound pressure level of left channel. Used for data
%                    display and analysis. Default value is 70.
%
%   See also: dietz2011interauralfunctions
%
%   References: dietz2011auditory

% AUTHOR: Mathias Dietz, Martin Klein-Hennig (implementation for AMT)

%   Copyright (C) 2002-2012   AG Medizinische Physik,
%                             Universitaet Oldenburg, Germany
%                             http://www.physik.uni-oldenburg.de/docs/medi
%
%   Authors: Tobias Peters (tobias@medi.physik.uni-oldenburg.de) 2002
%            Mathias Dietz (mathias.dietz@uni-oldenburg.de)      2006-2009
%            Martin Klein-Hennig (martin.klein.hennig@uni-oldenburg.de) 2011
 
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

% import default arguments from other functions
definput.import={'auditoryfilterbank','ihcenvelope'};
% changing default parameters
definput.importdefaults = { ...
    'flow',200, ...    % gammatone lowest frequency / Hz
    'fhigh',5000, ...  % gammatone highest frequency / Hz
    'basef',1000, ...  % auditory filter should be centered at basef / Hz
    'ihc_breebart' ... % use haircell parameters as in Breebarts model
};

% Gammatone filterbank parameters
%definput.keyvals.flow = 200;
%definput.keyvals.fhigh = 5000;
%definput.keyvals.basef = 1000;
%definput.keyvals.filters_per_ERB = 1;

% Preprocessing parameters
definput.keyvals.middle_ear_thr = [500 2000]; % Bandpass freqencies for middle ear transfer
definput.keyvals.middle_ear_order = 2;        % Only even numbers possible
definput.keyvals.compression_power = 0.4;
definput.keyvals.alpha = 0;                   % Internal noise strength
                                              % 65dB = 0.0354
% randn: add random noise with rms = alpha
% mini: set all values < alpha to alpha
definput.flags.int_noise_case = {'int_randn','int_mini'};

% === Binaural processor ===
% Parameters for filtering the haircell output
definput.keyvals.filter_order = 2;            % Used for both mod and fine
definput.keyvals.filter_attenuation_db = 10;  % Used for both mod and fine
% Finestructure filter
definput.keyvals.fine_filter_finesse = 3;
% Modulation filter
definput.keyvals.mod_center_frequency_hz = 135;
definput.keyvals.mod_filter_finesse = 8; % => bandwidth: 16.9 Hz

% For ild- or level-plugin
definput.keyvals.level_filter_cutoff_hz = 30;
definput.keyvals.level_filter_order = 2;

% Parameters for localization plugin
definput.keyvals.coh_param.max_abs_itd = 1e-3;
definput.keyvals.coh_param.tau_cycles  = 5;   % in cycles
definput.keyvals.coh_param.tau_s       = 10e-3; % in s for faller

% parameters for data display and analysis
definput.keyvals.signal_level_dB_SPL = 70; % sound pressure level of left channel

% display current process
displ                = 0;   % display current process (0 = do not)

[flags,kv]  = ltfatarghelper({},definput,varargin);



%% Model processing starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siglen=length(insig);
t = (0:siglen-1)/fs.';

%% ---- middle ear band pass filtering ------
% TODO: test middleear() function
if displ 
  disp(['band pass filtering of input according to middle ear transfer ' ...
        'charact. -> signal_me']); 
end
[b,a] = butter(kv.middle_ear_order,kv.middle_ear_thr(2)/(fs/2),'low');
inoutsig = filter(b,a,insig);
[b,a] = butter(kv.middle_ear_order,kv.middle_ear_thr(1)/(fs/2),'high');
inoutsig = filter(b,a,inoutsig);


%% ---- inner ear filterbank ------
if displ disp('splitting signal into frequency channels -> signal_filtered'); end
%[inoutsig,fc] = auditoryfilterbank(inoutsig,fs,'argimport',flags,kv);
% create filterbank
analyzer = gfb_analyzer_new(fs, kv.flow, kv.basef, kv.fhigh, kv.filters_per_ERB);
% center frequencies
fc = analyzer.center_frequencies_hz;
% get number of channels
channels = length(fc);

%% apply filterbank
inoutsig_left = gfb_analyzer_process(analyzer,inoutsig(:,1))';
inoutsig_right = gfb_analyzer_process(analyzer,inoutsig(:,2))';

% determine lowpass parameter
tau = kv.coh_param.tau_cycles./fc;

% cochlea compression to the power of 0.4
inoutsig_left = inoutsig_left.^kv.compression_power;
inoutsig_right = inoutsig_right.^kv.compression_power;

% rectification and lowpass filtering of filtered signals
% (haircell processing)
if displ disp('haircell processing of frequency bands -> hairc'); end
inoutsig_left = ihcenvelope(inoutsig_left,fs,'ihc_breebaart');
inoutsig_right = ihcenvelope(inoutsig_right,fs,'ihc_breebaart');
% only half-wave rectification for fine-structure channel
inoutsig_fine_left = max( inoutsig_left, 0 );
inoutsig_fine_right = max( inoutsig_right, 0 );

% adding internal noise
if flags.do_int_randn
  if displ disp('adding internal random noise -> hairc'); end
  addnoise = randn(size(inoutsig_left))*kv.alpha;
  inoutsig_left = inoutsig_left + addnoise;
  inoutsig_right = inoutsig_right + addnoise;
  inoutsig_fine_left = inoutsig_fine_left + addnoise;
  inoutsig_fine_right = inoutsig_fine_right + addnoise;
end
if flags.do_int_mini
  if displ disp('adding internal noise via minimum -> hairc'); end
  inoutsig_left = max(inoutsig_left,kv.alpha);
  inoutsig_right = max(inoutsig_right,kv.alpha);
  inoutsig_fine_left = max(inoutsig_fine_left,kv.alpha);
  inoutsig_fine_right = max(inoutsig_fine_right,kv.alpha);
end

%% === Binaural processor ===
%
% --- modulation filter ---
mod_filter_bandwidth_hz = kv.mod_center_frequency_hz/kv.mod_filter_finesse;
[inoutsig_mod_left, inoutsig_mod_right] = ... 
    gfb_envelope_filter(inoutsig_left, inoutsig_right, fs,...
    kv.mod_center_frequency_hz, mod_filter_bandwidth_hz, ...
    kv.filter_attenuation_db, kv.filter_order);
% binaural results for modulation filter
env = dietz2011interauralfunctions(...
    inoutsig_mod_left, inoutsig_mod_right,tau,kv.mod_center_frequency_hz+0*fc,...
    kv.signal_level_dB_SPL, kv.compression_power, kv.coh_param.tau_cycles, fs);
%
% --- fine structur filter ---
[inoutsig_fine_left, inoutsig_fine_right] =...
    gfb_envelope_filter(inoutsig_fine_left, inoutsig_fine_right, fs, fc,...
    fc/kv.fine_filter_finesse, kv.filter_attenuation_db, kv.filter_order);
% binaural results for fine structure filter 
fine = dietz2011interauralfunctions(...
    inoutsig_fine_left, inoutsig_fine_right, tau, fc,...
    kv.signal_level_dB_SPL, kv.compression_power, kv.coh_param.tau_cycles, fs);
% remove finestructure information > 1400 Hz
fine.f_inst(:,fc>1400)=[];
fine.ic(:,fc>1400)=[];
fine.ipd_lp(:,fc>1400)=[];
fine.itd_lp(:,fc>1400)=[];
%
% --- ILD ---
ild = ild_filter(inoutsig_left,inoutsig_right,kv.level_filter_cutoff_hz,...
                       kv.level_filter_order,kv.compression_power,fs);


if displ disp('finished'); end

end



%% gfb_envelope_filter %%%%%%%%%%%%%%%%%
function [envelopes_left, envelopes_right] = gfb_envelope_filter(insig_left, insig_right, sampling_rate_hz, center_frequency_hz,...
    bandwidth_hz, attenuation_db, gamma_filter_order)
% [envelopes_filtered, envelopes_sh_filtered] =...
%   gfb_envelope_filter(insig_left, insig_right, sampling_rate_hz, center_frequency_hz, bandwidth_hz, attenuation_db, gamma_filter_order);
%
% Filters each row of s1 and s2 with the gammatone filter defined by the input parameters.
% Takes both vectors and matrices.
%
% Input
%   insig_left, insig_right - signals to be filtered
%   sampling_rate_hz - sampling frequency / Hz
%   center_frequency_hz - centre frequency of the gammatone filter / Hz
%   bandwidth_hz - bandwidth of the gammatone filter at the level attenuation_db / Hz
%   attenuation_db - attenuation in dB at which the filter has bandwidth_hz
%   gamma_filter_order - order of the filter
    [N, M] = size(insig_left);
    if length(center_frequency_hz) == 1
        center_frequency_hz = center_frequency_hz * ones(1,M);
    end
    if isempty(bandwidth_hz) % default: width = 1 ERB
        recip_width1erb = diff(gfb_hz2erbscale(1:N/2));
        bandwidth_hz = round(1./recip_width1erb(round(center_frequency_hz)));
    elseif length(bandwidth_hz) == 1
        bandwidth_hz = bandwidth_hz * ones(1,M);
    end

    for ii=1:M
        analyzer = gfb_filter_new(sampling_rate_hz, center_frequency_hz(ii),...
            bandwidth_hz(ii), attenuation_db, gamma_filter_order);
        envelopes_left(:,ii)  = gfb_filter_process(analyzer, insig_left(:,ii)');
        envelopes_right(:,ii) = gfb_filter_process(analyzer, insig_right(:,ii)');
    end
end


%% ild_filter
function ild = ild_filter(insig_left,insig_right,lp_threshold_freq,lp_order,compression,fs)
    % calculation of ILD using the maxima of the envelopes
    % lowpass filtering
    [b,a] = butter(lp_order,lp_threshold_freq/(fs/2),'low');
    inoutsig_left  = filter(b,a,insig_left);
    inoutsig_right = filter(b,a,insig_right);
    ild = 20*log10(max(inoutsig_right,1e-4)./max(inoutsig_left,1e-4))/compression;
end
