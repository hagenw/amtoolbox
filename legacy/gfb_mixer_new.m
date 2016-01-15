function mixer = gfb_mixer_new(analyzer, delay, iterations)
%GFB_MIXER_NEW  Create new mixer
%   Usage: mixer = gfb_mixer_new(analyzer, delay)
%          mixer = gfb_mixer_new(analyzer, delay, iterations)
%
%   Input parameters:
%     analyzer : A gfb_analyzer structure as created by
%                gfb_analyzer_new. The mixer created by this function can
%                act as part of a synthesizer
%                that resynthesizes the output of this analyzer
%     delay :    A gfb_Delay structure as created by gfb_delay_new, together
%                with the mixer created by this function, this delay can
%                form a synthesizer that resynthesizes the output of the
%                analyzer
%     iterations : The gain factors are approximated numerically in
%                  iterations. If this parameter is omitted, then the
%                  number of iterations from analyzer will be used. 
%
%   `gfb_mixer_new` creates a `gfb_mixer` object with gain factors suitable
%   to calculate a weighted sum of the bands present in the output of the
%   given delay.  The gain factors are computed using a numerical optimization
%   method described in Herzke & Hohmann (2007).

% author   : tp 
% date     : Jan 2002, Nov 2003, Mar & Nov 2006, Jan Feb 2007

warning('Warning: GFB_MIXER_NEW will be removed in a future release. Use HOHMANN2002MIXER instead. ');

if (nargin < 4)
  iterations = analyzer.gaincalc_iterations;
end

mixer.type = 'gfb_mixer';
center_frequencies = analyzer.center_frequencies_hz;
number_of_bands = length(center_frequencies);
sampling_frequency = analyzer.fs;

% The center frequencies in the z plain
z_c = exp(2i * pi * center_frequencies(:) / sampling_frequency);

mixer.gains = ones(number_of_bands, 1);

% compute the frequency response of each filter (col) at the center
% frequencies of all filters (row)
  pos_f_response = hohmann2002freqz(analyzer, z_c);
  neg_f_response = hohmann2002freqz(analyzer, conj(z_c));

% apply delay and phase correction
for band = [1:number_of_bands]
  pos_f_response(:,band) = pos_f_response(:,band) * ...
    delay.phase_factors(band) .* ...
    z_c .^ -delay.delays_samples(band);
  neg_f_response(:,band) = neg_f_response(:,band) * ...
    delay.phase_factors(band) .* ...
    conj(z_c) .^ -delay.delays_samples(band);
end

% combine responses at positive and negative responses to yield
% responses for real part.
f_response = (pos_f_response + conj(neg_f_response)) / 2;

for i = 1:iterations
  % add selected spectrum of all bands with gain factors
  selected_spectrum = f_response * mixer.gains;

  % calculate better gain factors from result
  mixer.gains = mixer.gains ./ abs(selected_spectrum);
end
mixer.gains = mixer.gains.';
