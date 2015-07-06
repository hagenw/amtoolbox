% AMT - Single stages of auditory models
%
%  The AMT team, 2011 - 2014.
%
%  Peripheral stages
%     HEADPHONEFILTER    - FIR filter to model headphones 
%     MIDDLEEARFILTER    - FIR filter to model the middle ear
%     AUDITORYFILTERBANK - Linear auditory filterbank.
%     IHCENVELOPE        - Inner hair cell envelope extration.
%     ADAPTLOOP          - Adaptation loops.
%     KARJALAINEN1996    - Adaptation network.
%     DRNL               - Dual resonance non-linear filterbank.
%     MODFILTERBANK      - Modulation filter bank.
%     MODFILTERBANKEPSM  - Modulation filter bank from Joergensen model.
%
%  Binaural processing stages
%     LANGENDIJK2002COMP     - Comparison process from Langendijk 20002.
%     LINDEMANN1986BINCORR - Running cross-correlation between two signals.
%     EICELL             - Excitation-inhibition cell model by Breebaart.
%     CULLING2005BMLD    - BMLD calculation from Culling et al. (2005).
%     DIETZ2011UNWRAPITD - IPD to ITD transformation for the Dietz model
%     DIETZ2011FILTERBANK - filterbank of Dietz 2011 binaural model
%     DIETZ2011INTERAURALFUNCTIONS - Calculate interaural parameters for Dietz 2011 model
%     WIERSTORF2013ESTIMATEAZIMUTH - Azimuth position estimation based on dietz2011 or lindemann1986
%
%  The Takanen 2013 model   - Takanen 2013 binaural auditory model
%     TAKANEN2013CONTRACOMPARISON - Enhance contrast between hemispheres
%     TAKANEN2013CUECONSISTENCY - Check consistency before cue combination
%     TAKANEN2013DIRECTIONMAPPING - Map the directional cues to directions
%     TAKANEN2013FORMBINAURALACTIVITYMAP - Steer cues on a topographic map
%     TAKANEN2013LSO        - Model of the lateral superior olive
%     TAKANEN2013MSO        - Model of the medial superior olive
%     TAKANEN2013ONSETENHANCEMENT - Emphasize onsets on direction analysis
%     TAKANEN2013PERIPHERY  - Process input through the model of periphery
%     TAKANEN2013WBMSO      - Wideband medial superior olive model
%
%  Sagittal-plane localization models
%     BAUMGARTNER2013PMV2PPP - Calculate performance predictions from PMVs for baumgartner2013
%     BAUMGARTNER2014SPECTRALANALYSIS - Approximation of spectral analysis by auditory periphery
%     BAUMGARTNER2014GRADIENTEXTRACTION - Extraction of positive spectral gradients
%     BAUMGARTNER2014COMPARISONPROCESS - Comparison with direction-specific templates
%     BAUMGARTNER2014SIMILARITYESTIMATION - Similarity estimation with listener-specific sensitivity
%     BAUMGARTNER2014BINAURALWEIGHTING - Binaural combination of monaural similarity estimates
%     BAUMGARTNER2014SENSORIMOTORMAPPING - Response scatter induced by localization task
%     BAUMGARTNER2014CALIBRATION - Calibration of the model
%     BAUMGARTNER2014VIRTUALEXP - Performs a virtual sound-localization experiment
%     BAUMGARTNER2014PMV2PPP - Performance predictions from PMVs of baumgartner2014
%     BAUMGARTNER2014LIKELISTAT - Likelihood statistics for evaluation of model performance
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net

