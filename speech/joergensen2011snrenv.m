function [fcs_EPSM outSNRenvs ] = joergensen2011snrenv(MixEnv,NoiseEnv,fs)
%JOERGENSEN2011SNRENV  SNRenv of processed signal
%   Usage: [fcs_EPSM outSNRenvs ] = joergensen2011snrenv(MixEnv,NoiseEnv,fs);
%
%   Input parameters:
%     MixEnv       : The envelope of a signal which have been processed
%                    in some way
%     NoiseEnv     : The noise envelope with/without processing
%     fs           : Sampling frequency
%     outSNRenvs   : $1\times 7$ vector of overall SNRenv's in dB,
%                    one for each modulation filter center frequency
%     fcs_EPSM     : center-frequencies of the modulation filterbank
%     sEPSM_ExPtns : Envelope power excitation patterns for of the input
% 
%   `[fcs_EPSM,outSNRenvs]=joergensen2011snrenv(MixEnv,NoiseEnv,fs)`
%   calculates the SNRenv for the processed signal.

% AUTHOR: May 2012 Søren Jørgensen 

% Center frequencies of the modulation filters
% fcs = [1 2 4 8 16 32 64 ];

% Calculation of envelope power in 7 modulation filters
[fcs_EPSM Mix_ExcPtn]   =  modfilterbankepsm(MixEnv,0,fs);
[fcs_EPSM Noise_ExcPtn] =  modfilterbankepsm(NoiseEnv,0,fs);

% NaN values are set to zero
idx_nans =  isnan(Noise_ExcPtn);
for k = 1:length(idx_nans)
    if idx_nans(k) == 1
        Noise_ExcPtn = 0;
    end
end

% Noisefloor cannot exceed the mix, since they exist at the same time 
Noise_ExcPtn = min(Mix_ExcPtn,Noise_ExcPtn); 

% The noisefloor restricted to minimum 0.01 reflecting and internal noise
% threshold
Noise_ExcPtn = max(Noise_ExcPtn,0.01);
Mix_ExcPtn = max(Mix_ExcPtn,0.01);
% calculation of SNRenv
outSNRenvs = 10*log10(((Mix_ExcPtn-Noise_ExcPtn ) ./Noise_ExcPtn));

% SNRenv - values are truncated to minimum -30 dB.
outSNRenvs = max(-30,outSNRenvs);

% Excitation patterns
sEPSM_ExPtns = [Mix_ExcPtn'  Noise_ExcPtn'];


