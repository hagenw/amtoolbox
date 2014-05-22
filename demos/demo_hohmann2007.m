%DEMO_HOHMANN2007  Filterbank example
%
%   This example program demonstrates how to create and use an analysis
%   gammatone filterbank
%
% FIX ME: Not sure if this really implements Herzke & Hohmann (2007) or Hohmann (2002)

% author   : tp
% date     : Jan, Mar 2002, Nov 2006

flow = 70;
fhigh = 6700;
base_frequency_hz = 1000;
sampling_rate_hz = 16276;
filters_per_ERB = 1.0;

disp(['Building a filterbank for ', num2str(sampling_rate_hz), ...
      'Hz sampling frequency.']);
disp(['Lower cutoff frequency: ', num2str(flow), 'Hz']);
disp(['Upper cutoff frequency: ', num2str(fhigh), 'Hz']);
disp(['Base frequency        : ', num2str(base_frequency_hz), 'Hz']);
disp(['filters per ERB       : ', num2str(filters_per_ERB)]);
disp(' ')
analyzer = gfb_analyzer_new(sampling_rate_hz, flow, ...
                            base_frequency_hz, fhigh,...
			    filters_per_ERB);
bands = length(analyzer.center_frequencies_hz);
disp(['filterbank contains ', num2str(bands), ' filters:']);

fprintf(1,'%3s|%12s |%15s |%16s\n\n', ...
        '# ', 'f / Hz ', 'normalization', 'coefficient');

%% display filter parameters of the individual filters: %%
for band = 1:bands
  filter = analyzer.filters(band);
  fprintf(1,'%3d|%12f |%15e | %f + %fi\n', ...
	  band, analyzer.center_frequencies_hz(band), ...
	  filter.normalization_factor, ...
	  real(filter.coefficient), imag(filter.coefficient));
end


%%% plot the frequency response of the individual filters: %%%         

impulse = [1, zeros(1,8191)];                                          
[impulse_response, analyzer] = gfb_analyzer_process(analyzer, impulse);
frequency_response = fft(real(impulse_response)');                     
frequency = [0:8191] * sampling_rate_hz / 8192;                        

figure(1);
plot(frequency, 20 * log10(abs(frequency_response)));
axis([0, sampling_rate_hz/2, -40, 0]);
title('frequency response of the individual filters in this filterbank');
xlabel('frequency / Hz');
ylabel('filter response / dB');

disp(' ');
disp('Figure 1 shows the frequency response of the individual filters.');

