% AMT - Single stages of auditory models
%
%  The AMT team, 2011 - 2012.
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
%
%  Utility
%     PMV2PPP            - PMV to PPP conversion
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net
