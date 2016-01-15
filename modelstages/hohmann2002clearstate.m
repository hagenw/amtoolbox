function obj = hohmann2002clearstate(obj)
%HOHMANN2002  Reset states of a filter or a filterbank
%   Usage: filter = hohmann2002clearstate(filter)
%          fb = hohmann2002clearstate(fb)
%
%   `filter = hohmann2002clearstate(filter)` resets the states of the
%   filter created by |hohmann2002filter|
%
%   `fb = hohmann2002clearstate(fb)` resets the states of the
%   filterbank fb created by |hohmann2002|
% 
  
% Original author: Universitaet Oldenburg, tp (Jan 2002, Nov 2006, Feb 2007)
% Adapted to AMT: PM (Jan 2016)

if ~isfield(obj,'type'), error('Type of the object missing'); end
switch(obj.type)
  case 'gfb_Filter'
    obj.state = zeros(1, obj.gamma_order);
  case 'gfb_analyzer'
    for band = 1:length(obj.center_frequencies_hz)
        obj.filters(1, band).state = zeros(1, obj.filters(1, band).gamma_order);
    end
end
