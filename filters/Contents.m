% AMT - Filter functions
%
%  The AMT team, 2012 - 2014.
%
%  General routines
%     UFILTERBANKZ     - Apply multiple filters
%     FILTERBANKZ      - Apply multiple filters with non-equidistant downsampling
%     FILTERBANK_INIT  - Create control structure for FILTERBANK_BLOCK
%     FILTERBANK_BLOCK - Filterbank block processing
%
%  Auditory filters
%     GAMMATONE        - Gammatone filter coefficients
%     CQDFT            - FFT-based filter bank with constant relative bandwidth
%
%  Hohmann (2002) filterbank
%     GFB_ANALYZER_NEW            - Gammatone filterbank implementation, see demo_hohmann2012
%     GFB_ANALYZER_PROCESS        - All gfb_ functions are part of the hohmann2012 model
%     GFB_DELAY_NEW               - .
%     GFB_DELAY_PROCESS           - .
%     GFB_FILTER_NEW              - .
%     GFB_FILTER_PROCESS          - .
%     GFB_MIXER_NEW               - .
%     GFB_SYNTHESIZER_NEW         - .
%
%  Averaging
%     WEIGHTEDAVERAGEFILTER       - Part of the takanen2013 model
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net

