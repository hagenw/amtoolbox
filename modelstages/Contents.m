% AMT - Single stages of auditory models
%
%  The AMT team, 2011 - 2013.
%
%  Peripheral stages
%     HEADPHONEFILTER    - FIR filter to model headphones
%     MIDDLEEARFILTER    - FIR filter to model the middle ear
%     AUDITORYFILTERBANK - Linear auditory filterbank.
%     IHCENVELOPE        - Inner hair cell envelope extration.
%     ADAPTLOOP          - Adaptation loops.
%     DRNL               - Dual resonance non-linear filterbank.
%     MODFILTERBANK      - Modulation filter bank.
%
%  Binaural processing stages
%     LANGENDIJKCOMP     - Comparison process from Langendijk 20002.
%     LINDEMANNBINCORR   - Running cross-correlation between two signals.
%     EICELL             - Excitation-inhibition cell model by Breebaart.
%     CULLING2005BMLD    - BMLD calculation from Culling et al. (2005).
%     UNWRAP_ITD         - IPD to ITD transformation for the Dietz model
%     ITD2AZIMUTH        - ITD to azimuth conversion for Dietz and Lindemann model
%     LOOKUP_TABLE       - creating lookup tables for estimate_azimuth
%
%  The Takanen 2013 model   - XXX
%     TAKANEN2013CONTRACOMPARISON - XXX
%     TAKANEN2013CUECONSISTENCY - XXX
%     TAKANEN2013DIRECTIONMAPPING - XXX
%     TAKANEN2013FORMBINAURALACTIVITYMAP - XXX
%     TAKANEN2013LSO        - XXX
%     TAKANEN2013MSO        - XXX
%     TAKANEN2013ONSETENHANCEMENT - XX
%     TAKANEN2013PERIPHERY  - XXX
%     TAKANEN2013WBMSO      - XXX
%
%  Utility
%     PMV2PPP            - PMV to PPP conversion
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net
