function h = hohmann2002_freqz(obj, z)
%hohmann2002_freqz  Frequency response of hohmann2002 filter or filterbank
%   Usage: h = hohmann2002_freqz(filter, z)
%          H = hohmann2002_freqz(fb, z)
%
%   `h = hohmann2002_freqz(filter, z)` returns the frequency response of
%   filter created by |hohmann2002_filter|.
%    
%   `h = hohmann2002_freqz(fb, z)` returns the frequency responses of the
%   individual filters in the filterbank fb created by |hohmann2002|. The
%   responses of the individual filters are stored in columns of h. 
%
%   Input parameters:
%     filter  : A filter created by |hohmann2002_filter|.
%     fb      : A filterbank created by |hohmann2002|.
%     z       : A vector of frequencies in z-plane for which the frequency 
%               response will be computed. $z = \exp(2i\cdot pi\cdot f/fs)$
%
%   Output parameters:
%     h       : The complex frequency response at z. Each column represents
%               a response of a filter.
%
%   See also: exp_hohmann2002 demo_hohmann2002
%
%   References: hohmann2002
%

% author   : Universitaet Oldenburg, tp (Jan & Nov 2006, Jan Feb 2007)
% Adapted to AMT (PM, Jan 2016) from functions gfb_*_zresponse

if ~isfield(obj,'type'), error('Type of the object missing'); end
switch(obj.type)
  case 'gfb_Filter'
    h = (1 - obj.coefficient ./ z) .^ -obj.gamma_order * obj.normalization_factor;
  case 'gfb_analyzer'
    number_of_bands = length(obj.center_frequencies_hz);
    z = z(:);
    h = ones(length(z), number_of_bands);

    for band = 1:number_of_bands
      filter = obj.filters(band);
      h(:,band) = (1 - filter.coefficient ./ z) .^ -filter.gamma_order * filter.normalization_factor;
    end
  otherwise
    error('Unknown type of HOHMANN2002 filter object');
end

