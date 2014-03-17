Authors: Michael Weinert (Michael.Weinert@rub.de), Gerald Enzner (Gerald.Enzner@rub.de)
Date: 	 21-01-2012


These tools will demonstrate signal processing for continuous-azimuth head-related impulse response (HRIR) acquisition [Enzner2008]. Based on a given set of example recordings, you are able to create a dense HRIR data set. In this example, the horizontal plane (elevation = 0 °) was considered. The example measurements used two different measurement stimuli: white noise and perfect sweeps [Antweiler2012]. The subject was a HEAD acoustics dummy head KMS 3.II. In order to make use of the computed HRIRs, the HRIRs can be spatially sampled at arbitrary azimuth-resolution. A simple binaural rendering (demo_enzner2008.m) will demonstrate immediate usage of the computed HRIRs.   

The perfect sweep signal used in this demonstration can be generated with a separate tool (perfectsweep.m).

The following files are provided with the tool-chain:

List of m-files (tools):
------------------------

\hrtf\continuous-azimuth HRIR\
	enzner2008.m	  : Computes HRIRs with continuous-azimuth HRIR representation and stores data at abitrary azimuth-resolution

\experiments\
	exp_enzner2008.m  : Creates figures like [Enzner2008, Fig. 2], [Enzner2009, Fig. 4]

\demos\demo_enzner2008\
	demo_enzner2008.m : Binaural rendering to demonstrate the usage of the computed data


Perfect sweep generator:
------------------------

\signals\PerfectSweep\
	
	perfectsweep.m	: Generates a perfect sweep sequence (one period of the periodic perfect sweep, normalized to the range -1...1)


Example files:
--------------

1. Measurement stimulus: white noise

   loudspeaker driving signal:			\hrtf\continuous-azimuth HRIR\signals\playback\example_1ch_white_noise_playback.wav   
   recorded earsignals (stero):	   		\hrtf\continuous-azimuth HRIR\measurement\example_1ch_white_noise_earsignals.wav
   reference (mono recording without head):  	\hrtf\continuous-azimuth HRIR\signals\reference\example_1ch_white_noise_reference.wav


2. Measurement stimulus: perfect sweeps (N=308)

   loudspeaker driving signal:			\hrtf\continuous-azimuth HRIR\signals\playback\example_1ch_PSWEEP308_playback.wav   
   recorded earsignals (stero):			\hrtf\continuous-azimuth HRIR\measurement\example_1ch_PSWEEP308_earsignals.wav
   reference (mono recording without head):	\hrtf\continuous-azimuth HRIR\signals\reference\example_1ch_PSWEEP308_reference.wav


3. Input signal example for rendering demo (speech):

   \hrtf\continuous-azimuth HRIR\EBU_SQAM49.wav (see http://tech.ebu.ch/publications/sqamcd)


Additional information about NLMS Input signals
-----------------------------------------------


A) Using loudspeaker driving signals:

The easyest way to compute impulse responses using the NLMS-algorithm is to use the algorithm with the loudspeaker driving signals as input.

Keep in mind, that the computed impulse responses then also contain the characteristics of the loudspeaker and the distance from the lousdpeaker to the in-ear microphones. To take this into consideration,  use sufficiently large filter lengths "h_length" or choose a negative value for "sys_latency".


B) Using reference signals measured at the head position without head:

Reference signals are recorded loudspeaker driving signals at the position of the listeners head, while the listener is absent. The general measurement setup has to be the same as used when performing the measurement of a subject, e.g., loudness, calibration of amplifieres, etc. Thus, reference signals contain some characteristics of the measurement system, e.g., loustpeaker and amplifiers. Using reference signals, those characteristics can be eliminated from the HRIR results during the computation. 

Hence, the given reference recordings can only be used in combination with the given examples of measurements (recorded ear signals). To use reference signals in your own measurements, you have to create your own reference signals. The given loudspeaker driving signals can be used for this purpose.

Keep in mind, that you have to ensure causality of HRIRs by adjusting the parameter "sys_latency" to a suitable value (see file RUBHRIR_NLMS.m)!  

