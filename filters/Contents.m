% AMT - Filter functions
%
%  The AMT team, 2012.
%
%  General routines
%     UFILTERBANKZ     - Apply multiple filters
%     FILTERBANKZ      - Apply multiple filters with non-equidistant downsampling
%     FILTERBANK_INIT  - Create control structure for FILTERBANK_BLOCK
%     FILTERBANK_BLOCK - Filterbank block processing
%     AVERAGINGFB      - Averaging rectangular filter bank
%
%  Auditory filters
%     GAMMATONE        - Gammatone filter coefficients
%
%  Hohmann (2002) filterbank
%     GFB_ANALYZER_CLEAR_STATE    - XXX 
%     GFB_ANALYZER_NEW            - XXX 
%     GFB_ANALYZER_PROCESS        - XXX 
%     GFB_ANALYZER_ZRESPONSE      - XXX 
%     GFB_DELAY_CLEAR_STATE       - XXX 
%     GFB_DELAY_NEW               - XXX 
%     GFB_DELAY_PROCESS           - XXX 
%     GFB_FILTER_CLEAR_STATE      - XXX 
%     GFB_FILTER_NEW              - XXX 
%     GFB_FILTER_PROCESS          - XXX 
%     GFB_FILTER_ZRESPONSE        - XXX 
%     GFB_MIXER_NEW               - XXX 
%     GFB_MIXER_PROCESS           - XXX 
%     GFB_SET_CONSTANTS           - XXX 
%     GFB_SYNTHESIZER_CLEAR_STATE - XXX 
%     GFB_SYNTHESIZER_NEW         - XXX 
%     GFB_SYNTHESIZER_PROCESS     - XXX 
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net

