%DEMO_HOHMANN2002  Filterbank example
%
%   This example program demonstrates how to create and use the
%   gammatone filterbank within the framework of Hohmann (2002).
%
%   See also: exp_hohmann2002 hohmann2002 hohmann2002process

% author   : tp
% date     : Jan, Mar 2002, Nov 2006

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

%% display filter parameters of the individual filters
for band = 1:bands
  filter = analyzer.filters(band);
  fprintf(1,'%3d|%12f |%15e | %f + %fi\n', ...
	  band, analyzer.center_frequencies_hz(band), ...
	  filter.normalization_factor, ...
	  real(filter.coefficient), imag(filter.coefficient));
end


%% plot the frequency response of the individual filters     

frequency = 0:2:fs/2;
h=hohmann2002freqz(analyzer,exp(2*1i*pi*frequency/fs));

figure(1);
plot(frequency, 20 * log10(abs(h)));
axis([0, fs/2, -40, 7]);
title('frequency response of the individual filters in this filterbank');
xlabel('frequency / Hz');
ylabel('filter response / dB');

amtdisp(' ');
amtdisp('Figure 1 shows the frequency response of the individual filters.');

