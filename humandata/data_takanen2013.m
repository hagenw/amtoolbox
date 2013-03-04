%DATA_TAKANEN2013 Data applied in the model by Takanen, Santala and Pulkki
%
%   XXX This should be a function and not a script
%
%   This script describes and loads the pre-computed data structures that
%   are applied in the different processing steps of the count-comparison
%   based binaural auditory model by Takanen, Santala and Pulkki.
%
%   See also: takanen2013
%
%   References: takanen2013a

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland

%% Data employed in the model of periphery
%frequency-dependent delays of the cochlea model by Verhulst in the
%characteristic frequencies of the model by Takanen et al. are employed
%in the *takanen2013periphery.m*-function. The delays were determined by
%analyzing the the cross-correlation-functions between different
%frequency bands.
load takanen2013_cochleardelays.mat -mat

%% Data employed in the MSO model
%in the *takanen2013mso.m*-function the ipsilateral input signals are
%divided by scaling values and thereafter limited between 0 and 1. The
%scaling values were obtained by computing the average level of the
%*periphOutput.left* for a pink noise signal reproduced at 30 dB SPL
load takanen2013_msolimits.mat -mat

%% Data employed in the Wide-band MSO model
%in the *takanen2013widebandmso.m*-function the energy output is multiplied
%so that the energy of the wide-band MSO model is in a similar level as
%compared to the energies of the narrowband MSO and LSO models. The values
%with which the energy output is multiplied were obtained through an
%iterative process.
load takanen2013_wbmsomultp.mat -mat

%% Data employed in the direction mapping
%in the *takanen2013directionmapping.m*, the outputs of the different MSO
%and LSO models are mapped into directions ranging from -20 to 90. The
%mapping is implemented following the idea of self-calibration, where the
%outputs of the models are compared separately to sets of reference values
%computed with HRTF-processed samples. The sample employed for the
%narrowband MSO and the LSO models was a 80-ms-long pink noise burst,
%whereas the sample for the wide-band MSO model was a impulse response of a
%first-order Butterworth lowpass filter with a cut-off frequency of 500 Hz
load takanen2013_lookuptable.mat -mat

%% Data employed in the onset contrast enhancement
%in the *takanen2013onsetenhancement.m*-function, the energies of the
%short-term and of the long-term directional cues are scaled to a similar
%level with the help of a set of pre-computed values, the values which were
%obtained by the computing the average levels of the two energies for a
%pink noise burst convolved with binaural room impulse response of a
%concert hall.
load takanen2013_onsetmultp.mat -mat

%% Data employed in the forming of the binaural activity map
%in the *takanen2013formbinauralactivitymap.m*-function the levels of the
%what cues are scaled in order to visualize the evoked activations on the 
%binaural activity map. The values were obtained by computing the average
%levels of the what cue for a pink noise burst reproduced at 60 dB SPL.
load takanen2013_periphenergyaverages.mat -mat