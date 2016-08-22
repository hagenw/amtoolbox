% AMT - Single stages of auditory models
%
%  The AMT team, 2011 - 2014.
%
%  HRTF processing
%     ZIEGELWANGER2014ONAXIS      - On-axis version of Ziegelwanger et al. (2014)
%     ZIEGELWANGER2014OFFAXIS     - Off-axis version of Ziegelwanger et al. (2014)
%     ZIEGELWANGER2013ONAXIS      - On-axis version of Ziegelwanger et al. (2013) (legacy only)
%     ZIEGELWANGER2013OFFAXIS     - Off-axis version of Ziegelwanger et al. (2013) (legacy only)
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
%     ZIEGELWANGER2014OFFAXIS - synthesis of TOAs for the off-axis model
%     ZIEGELWANGER2014ONAXIS  - synthesis of TOAs for the on-axis model
%     LANGENDIJK2002LIKELIHOOD - Likelihood estimation
%
%  Preprocessing stages
%     BREEBAART2001PREPROC  - EI-cell output
%
%     LINDEMANN1986CENTROID - Centroid of the cross-correlation activation
%
%     ROENNE2012CHIRP          - Simulate chirp evoked ABRs
%     ROENNE2012CLICK          - Simulate ABR respone to click
%     ROENNE2012TONEBURSTS     - Simulate tone burst evoked ABR wave V latencies
%
%   Gammatone filterbank framework from Hohmann (2002)
%     HOHMANN2002DELAY            - Create a delay object
%     HOHMANN2002MIXER            - Create a mixer object
%     HOHMANN2002SYNTH            - Create a synthesis object 
%     HOHMANN2002FILT             - Create a single Gammatone filter object
%     HOHMANN2002PROCESS          - Process the input signals by the corresponding filterbank object
%     HOHMANN2002CLEARSTATE       - Clears the state of the filterbank object
%     HOHMANN2002FREQZ            - Calculates frequency response of a filterbank object
%     HOHMANN2002CLEARSTATE       - Clears the state of the filterbank objects
%     HOHMANN2002FREQZ            - Calculates frequency response of a filterbank object
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
%  Binaural masking level differences from Breebaart et al. (2001)
%     BREEBAART2001PREPROC  - Preprocessor of a binaural signal with EI-cell output
%     BREEBAART2001CENTRALPROC - Central processor taking decision in an experiment
%     BREEBAART2001SIGGEN - Generates signals required for an experiment
% 
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

