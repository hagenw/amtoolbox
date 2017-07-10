% AMT - General functions used within the AMT
%
%  Signal levels and thresholds
%     DBSPL            - SPL of signal measured in dB.
%     SETDBSPL         - Specify SPL of signal.
%     ABSOLUTETHRESHOLD    - Absolute threshold of hearing
%
%  Emulation of experiments
%     EMUEXP         - Emulate psychoacoustic experiments
%
%  Plotting
%     AUDSPECGRAM      - Auditory spectrogram.
%     MODSPECGRAM      - Temporal modulation spectrogram
%     STMODSPECGRAM    - Spectro-temporal modulation spectrogram
%
%  Filters
%     GAMMATONE        - Calculate Gammatone filter coefficients
%     AUDITORYFILTERBANK - Linear auditory filterbank.
%     IHCENVELOPE        - Inner hair cell envelope extration.
%     ADAPTLOOP          - Adaptation loops.
%     MODFILTERBANK      - Modulation filter bank.
%     HEADPHONEFILTER    - FIR filter to model headphones 
%     MIDDLEEARFILTER    - FIR filter to model the middle ear
%     UFILTERBANKZ     - Apply multiple filters
%     FILTERBANKZ      - Apply multiple filters with non-equidistant downsampling
%
%  Sound localization
%     LOCALIZATIONERROR - Calculates various localization errors from localization responses
%     ITD2ANGLE        - Convert ITD to an angle using a lookup table
%     itd2angle_lookuptable - Create the lookup table
%
%  Speech
%     SIIWEIGHTINGS        - Speech intelligibility weighted by frequency