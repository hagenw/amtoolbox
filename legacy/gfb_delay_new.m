function delay = gfb_delay_new(analyzer, delay_samples)
%GFB_DELAY_NEW  Create new delay object
%   Usage: delay=gfb_delay_new(analyzer,delay_samples)
%
%   Input parameters:
%     analyzer      : The gfb_analyzer structure as returned by
%                     |gfb_analyzer_new|.
%     delay_samples : The desired group delay in samples. Must be at least 1,
%                     because of the way the phase factors are computed. Larger
%                     delays lead to better signal quality
%   Output parameters:
%     delay         : The new gfb_delay object
%
%   `gfb_delay_new(analyzer, delay_samples)` creates a new `gfb_delay
%   object` that can act as the first stage of a synthesizer that
%   resynthesizes the output of the gammatone filterbank analyzer.  The
%   purpose of the delay object is to delay the output of each band by a
%   band-dependent ammount of samples, so that the envelope of the impulse
%   response of the analyzer is as large as possible at the desired delay.
%   Additionally, the delay object will multiply this delayed output with a
%   band-dependent complex phase factor, so that the real part of the
%   impulse response has a local maximum at the desired delay.  Finally, the
%   delay object will output only the real part of each band.
%  
%   The phase factors are approximated numerically in this constructor,
%   using a method described in Herzke & Hohmann (2007).  The
%   approximation assumes parabolic behaviour of the real part of the
%   impulse response in the region of the desired local maximum: The phase
%   factors are chosen so that the real parts of the impulse response in
%   the samples directly preceeding and following the desired local
%   maximum will be equal after multiplication with the pase factor.
%
%   References: herzke2007improved

% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002; Nov 2003; Mar Jun Nov 2006; Jan Feb 2007

warning('Warning: GFB_DELAY_NEW will be removed in a future release. Use hohmann2002_delay instead. ');

delay.type           = 'gfb_Delay';

  analyzer             = hohmann2002_clearstate(analyzer);
  impulse              = zeros(1, delay_samples + 2);
  impulse(1)           = 1;

    impulse_response = ...
      gfb_analyzer_process(analyzer, impulse);

number_of_bands      = size(impulse_response, 1);

[dummy, max_indices] = max(abs(impulse_response(:,1:(delay_samples+1))).');

delay.delays_samples = delay_samples + 1 - max_indices;

delay.memory         = zeros(number_of_bands, max(delay.delays_samples));

slopes = zeros(1, number_of_bands);
for band = [1:number_of_bands]
  band_max_index = max_indices(band);
  slopes(band) = (impulse_response(band, band_max_index+1) - ...
                  impulse_response(band, band_max_index-1));
end
slopes = slopes ./ abs(slopes);
delay.phase_factors = 1i ./ slopes;

