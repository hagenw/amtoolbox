% AMT - General functions
%
%  The AMT team, 2012 - 2014.
%
%  General plotting
%     AUDSPECGRAM      - Auditory spectrogram.
%     MODSPECGRAM      - Temporal modulation spectrogram
%     STMODSPECGRAM    - Spectro-temporal modulation spectrogram
%
%   Auditory filters
%     GAMMATONE        - Calculate Gammatone filter coefficients
%     CQDFT            - FFT-based filter bank with constant relative bandwidth
%     AUDITORYFILTERBANK - Linear auditory filterbank.
%     IHCENVELOPE        - Inner hair cell envelope extration.
%     ADAPTLOOP          - Adaptation loops.
%     DRNL               - Dual resonance non-linear filterbank.
%     MODFILTERBANK      - Modulation filter bank.
%
%   Averaging
%     WEIGHTEDAVERAGEFILTER       - Part of the takanen2013 model
%
%  ITD
%     ITD2ANGLE        - Convert ITD to an angle using a lookup table
%     ITD2ANGLELOOKUPTABLE - Create the lookup table
%
%  Signal levels
%     DBSPL            - SPL of signal measured in dB.
%     SETDBSPL         - Specify SPL of signal.
%
%  Tools
%     LOCALIZATIONERROR - Calculates various localization errors from localization responses
%     HEADPHONEFILTER    - FIR filter to model headphones 
%     MIDDLEEARFILTER    - FIR filter to model the middle ear
%     MODFILTERBANKEPSM  - Modulation filter bank (model from Joergensen et al., 2012).
%     EICELL             - Excitation-inhibition cell (model by Breebaart et al., 2001).
%
%   General routines
%     UFILTERBANKZ     - Apply multiple filters
%     FILTERBANKZ      - Apply multiple filters with non-equidistant downsampling

%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net

