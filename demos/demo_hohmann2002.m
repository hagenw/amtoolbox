%DEMO_HOHMANN2002  Shows how to use the gammatone filterbank from Hohmann(2002)
%
%   Part I: This example creates a 4th order gammatone filter with a center
%   frequency of 1000Hz and a 3dB-bandwidth of 100Hz, suitable for 
%   input signals with a sampling frequency of 10kHz.
%
%   Part II: This example program demonstrates how to create and use the
%   gammatone filterbank within the framework of Hohmann (2002).
%
%   Part III: This Example demonstrates how to create and how to use the
%   combined analysis-synthesis Filterbank system.
%
%   See also: exp_hohmann2002 hohmann2002 hohmann2002process

% author   : tp
% date     : Jan, Mar 2002, Nov 2006
% last modified: Robert Baumgartner, Jan 18, 2016

%% Part I: Filter example

%%% create the example filter: %%%

center_frequency_hz   =  1000;
bandwidth_hz          =   100;
attenuation_db        =     3;
fs                    = 44100;
filter_order          =     4;

filter = hohmann2002filter(fs, center_frequency_hz, ...
                        bandwidth_hz, attenuation_db, filter_order);

%%% print the filter's parameters to the screen %%%

amtdisp(' ');
amtdisp(['The filter coefficient of this filter is: ', ...
      num2str(real(filter.coefficient)),            ...
      ' + ',                                        ...
      num2str(imag(filter.coefficient))]);
amtdisp(['Its normalization factor is             : ', ...
      num2str(filter.normalization_factor)]);


%%% plot the impulse response and the frequency response of this filter: %%%
impulse_samples            = 8192;
impulse_response_samples   =  800;
impulse                    = [1, zeros(1,impulse_samples - 1)];

[impulse_response, filter] = hohmann2002process(filter, impulse);

figure(1);
plot([0:impulse_response_samples-1], ...
         [real(impulse_response(1:impulse_response_samples)); ...
          imag(impulse_response(1:impulse_response_samples)); ...
          abs(impulse_response(1:impulse_response_samples))]);
axis([0,impulse_response_samples, -0.035,0.035]);
title('impulse response of example gammatone filter');
xlabel('sample number');
ylabel('filter output');

amtdisp(' ');
amtdisp(['Figure 1 shows the first ',num2str(impulse_response_samples),' samples of']);
amtdisp('the impulse response of a 4th order gammatone filter with a center'); 
amtdisp(['frequency of ',center_frequency_hz,'Hz and a 3dB-bandwidth of ',bandwidth_hz,'Hz.']);
amtdisp('Real part, imaginary part, and absolute value of the impulse ');
amtdisp('response are plotted as lines 1, 2, and 3, respectively.     ');

%%% plot the frequency response of this filter: %%%

frequency_response = abs(fft(real(impulse_response.')));
figure(2);
plot([0:(impulse_samples - 1)] *                     ...
     (fs / impulse_samples),       ...
     20 * log10(frequency_response));
axis([0,1500, -40,0]);
title('frequency response of example gammatone filter');
xlabel('frequency / Hz');
ylabel('filter response / dB');

amtdisp(' ');
amtdisp('Figure 2 shows the frequency response function of this filter');
amtdisp('in dB over frequency in Hz.');

%% Part II: Filterbank example

flow = 70;
fhigh = 6700;
base_frequency_hz = 1000;
fs = 16276;
filters_per_ERB = 1.0;

amtdisp(['Building a filterbank for ', num2str(fs), ...
      'Hz sampling frequency.']);
amtdisp(['Lower cutoff frequency: ', num2str(flow), 'Hz']);
amtdisp(['Upper cutoff frequency: ', num2str(fhigh), 'Hz']);
amtdisp(['Base frequency        : ', num2str(base_frequency_hz), 'Hz']);
amtdisp(['filters per ERB       : ', num2str(filters_per_ERB)]);
amtdisp(' ')
analyzer = hohmann2002(fs, flow, base_frequency_hz, fhigh, filters_per_ERB);
bands = length(analyzer.center_frequencies_hz);
amtdisp(['filterbank contains ', num2str(bands), ' filters:']);

fprintf(1,'%3s|%12s |%15s |%16s\n\n', ...
        '# ', 'f / Hz ', 'normalization', 'coefficient');

%%% display filter parameters of the individual filters
for band = 1:bands
  filter = analyzer.filters(band);
  fprintf(1,'%3d|%12f |%15e | %f + %fi\n', ...
	  band, analyzer.center_frequencies_hz(band), ...
	  filter.normalization_factor, ...
	  real(filter.coefficient), imag(filter.coefficient));
end


%%% plot the frequency response of the individual filters     

frequency = 0:2:fs/2;
h=hohmann2002freqz(analyzer,exp(2*1i*pi*frequency/fs));

figure(3);
plot(frequency, 20 * log10(abs(h)));
axis([0, fs/2, -40, 7]);
title('frequency response of the individual filters in this filterbank');
xlabel('frequency / Hz');
ylabel('filter response / dB');

amtdisp(' ');
amtdisp('Figure 3 shows the frequency response of the individual filters.');


%% Part III: Example for combined analysis-synthesis filterbank system

%%% First, create a filterbank analyzer as in Part I %%%

flow =    70;
fhigh =  6700;
base_frequency_hz         =  1000;
sampling_rate_hz          = 16276;
filters_per_ERB           =     1.0;
desired_delay_in_seconds  =     0.004;
filter_order              =     4;
bandwidth_factor          =     1.0;

amtdisp('Building analysis filterbank');
analyzer = hohmann2002(sampling_rate_hz, flow, ...
                            base_frequency_hz, fhigh,...
			    filters_per_ERB, filter_order, bandwidth_factor);


%%% Now create a synthesizer that can resynthesize the analyzer's output %%%

amtdisp(['Building synthesizer for an analysis-synthesis delay of ', ...
      num2str(desired_delay_in_seconds), ' seconds']);
synthesizer = hohmann2002synth(analyzer, desired_delay_in_seconds);

%%% Extract the synthesizer's parameters %%%
amtdisp(['The synthesizers parameters:';...
      '----------------------------']);
delay = synthesizer.delay;
mixer = synthesizer.mixer;

bands = length(mixer.gains);

fprintf(1,'%3s|%7s | %22s | %5s\n\n', ...
        '# ', 'delay ', 'phase factor    ', 'gain / dB');

for band = 1:bands
  fprintf(1,'%3d|%7d | %9f + %9fi | %5.2f\n', ...
	  band,                            ...
	  delay.delays_samples(band),      ...
	  real(delay.phase_factors(band)), ...
	  imag(delay.phase_factors(band)), ...
	  20*log10(mixer.gains(band)));
end

%%%  plot the resynthesized impulse and the frequency response of the  %%%
%%%  analysis-synthesis system                                         %%%

impulse = [1, zeros(1,8191)];                                          
[analyzed_impulse, analyzer] = hohmann2002process(analyzer, impulse);
[resynthesized_impulse, synthesizer] = ...
    hohmann2002process(synthesizer, analyzed_impulse);

figure(4);
plot([0:8191]/sampling_rate_hz*1e3, resynthesized_impulse);
axis([40/sampling_rate_hz*1e3, 120/sampling_rate_hz*1e3, -1, 1]);
title('impulse response of the analysis-synthesis system');
xlabel('time / ms');
ylabel('system output');

amtdisp(' ');
amtdisp('Figure 4 shows the impulse response of the analysis-synthesis');
amtdisp('system in the time domain.');

frequency = [0:8191] * sampling_rate_hz / 8192;
figure(5)
plot(frequency, 20 * log10(abs(fft(resynthesized_impulse'))));
axis([0, sampling_rate_hz/2, -40, 5]);
title('frequency response of the analysis-synthesis-system');
xlabel('frequency / Hz');
ylabel('system response level / dB'); 

amtdisp(' ');
amtdisp('Figure 5 shows its frequency response.');