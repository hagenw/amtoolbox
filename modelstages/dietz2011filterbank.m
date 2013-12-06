function [outsig_fine,fc_fine,outsig_env,fc_env,outsig_ild] = dietz2011filterbank(insig,fs,fc,varargin)
%DIETZ2011FILTERBANK  filterbank of Dietz 2011 binaural model
%   Usage: [...] = dietz2011filterbank(insig,fs,fc);
%
%   Input parameters:
%       insig       : binaural signal for which values should be calculated
%       fs          : sampling rate (Hz)
%       fc          : center frequencies of gammatone filterbank
%
%   Output parameters:
%       outsig_fine : output signal of fine structure filter
%       fc_fine     : center frequencies processed with fine structure filter
%       outsig_env  : output signal of envelope filter
%       fc_env      : center frequencies processed with modulation filter
%       outsig_ild  : output signal of ILD filter
%
%   `dietz2011filterbank(insig,fs,fc)` filters all frequency channels
%   of insig with a modulation filterbank consisting of three filters.
%   One centered at the center frequencies for the fine structure of the
%   signals. One centered at a fixed frequency of 135 Hz for the envelope of
%   the signals. And one a lowpass filter with a cutoff frequency of 30 Hz
%   for the calculation of the interaural level difference.
%
%   `dietz2011filterbank` accepts the following optional parameters:
%
%     'filter_order',fo
%                    Filter order for the two gammatone filter used for the fine
%                    structure and envelope of the modulation filter bank. The
%                    default value is 2.
%
%     'filter_attenuation_db',fadb
%                    Filter attenuation for the two gammatone filter used for the fine
%                    structure and envelope of the modulation filter bank. The
%                    default value is 10.
%                     
%     'fine_filter_finesse',fff
%                    Filter finesse (determines the bandwidth with fc/finesse)
%                    for the fine structure gammatone filter. The defulat value
%                    is 3.
%
%     'mod_center_frequency_hz',mcf_hz
%                    Center frequency of the gammatone envelope filter. The
%                    default value is 135.
%
%     'mod_filter_finesse',mff
%                    Filter finesse (determines the bandwidth with fc/finesse)
%                    for the envelope gammatone filter. The defulat value is 8.
% 
%     'level_filter_cutoff_hz',lfc_hz
%                    Cutoff frequency off the low pass filter used for ILD
%                    calculation. The default value is 30.
%
%     'level_filter_order',lforder
%                    Order of low pass filter for the ILD calculation. The
%                    default value is 2.
%
%   See also: dietz2011, dietz2011interauralfunctions
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
 
  
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig)
  error('%s: insig has to be a numeric signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

if ~isnumeric(fc)
  error('%s: fc has to be a numeric signal!',upper(mfilename));
end

definput.import = {'dietz2011filterbank'};
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Model processing starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Split the input signal:
% - below 1.4 kHz the fine structure filter is applied
% - above 1.4 kHz the modulation filter is applied
fc_fine = fc(fc<=1400);
fc_env = fc(fc>1400);
insig_fine = insig(:,fc<=1400,:);
insig_env = insig(:,fc>1400,:);

% --- fine structur filter ---
outsig_fine = zeros(size(insig_fine));
% gammatone filter centered at the center frequency for every frequency channel
for ii=1:length(fc_fine)
  gammatone_filter = gfb_filter_new(fs, fc_fine(ii), ...
    fc_fine(ii)/kv.fine_filter_finesse, ...
    kv.filter_attenuation_db, kv.filter_order);
  outsig_fine(:,ii,1) = gfb_filter_process(gammatone_filter, insig_fine(:,ii,1)');
  outsig_fine(:,ii,2) = gfb_filter_process(gammatone_filter, insig_fine(:,ii,2)');
end

% --- modulation/envelope filter ---
outsig_env = zeros(size(insig_env));
% gammatone filter centered at a fixed frequency for every frequency channel
gammatone_filter = gfb_filter_new(fs, kv.mod_center_frequency_hz, ...
  kv.mod_center_frequency_hz/kv.mod_filter_finesse, ...
  kv.filter_attenuation_db, kv.filter_order);
for ii=1:length(fc_env)
  outsig_env(:,ii,1) = gfb_filter_process(gammatone_filter, insig_env(:,ii,1)');
  outsig_env(:,ii,2) = gfb_filter_process(gammatone_filter, insig_env(:,ii,2)');
end

% --- ILD filter ---
% low pass filter with a fixed cutoff frequency for every frequency channel
[b,a] = butter(kv.level_filter_order,kv.level_filter_cutoff_hz/(fs/2),'low');
outsig_ild = filter(b,a,insig);

% vim: set sw=2 ts=2 expandtab textwidth=80: 
