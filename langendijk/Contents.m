% AMT - Langendijk's localization model
%
%  The AMT team, 2011.
%
%  Model
%     langendijk       - Langendijk & Bronkhorst 2002 auditory model 
%  Related routines
%     plotlangendijk   - Plot the pdf-matrices with grey colormap
%     likelilangendijk - Likelihood estimation to quantify model performance
%     plotlikelilangendijk - Plot the likelihood statistics 
%
%  General routines
%     langecomp        - Comparison process used in Langendijk's model
%     discreteinvrnd   - Invert a discrete distribution with probability p
%     halfconv         - Fast convolution, one way (fft but without ifft)
%     bmp2gr           - Convert a bitmap to gain response data (in dB)
%     gr2ir            - Convert gain responses to impulse responses
%
%  Routines and provided data, necessary for demo_langendijk
%     langendijk2002-P3 - DTF data and response patterns of listener P3
%     langendijk2002-P6 - DTF data and response patterns of listener P6
%     data_langendijk2002 - Returns response patterns and DTF data of
%                           given figures in Langendijk's paper
%     hrtf_M_langendijk2002 P3 - DTF data of P3 converted to ARI format  
%     hrtf_M_langendijk2002 P6 - DTF data of P6 converted to ARI format 
%     LangendijkColormap - Special colormap to get comparable results to
%                          the results in the paper of Langendijk & Bronkhorst 2002
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net