function zresponse = gfb_filter_zresponse(filter, z)
%GFB_FILTER_ZRESPONSE  Filter response at freqenzy z
%   Usage: zresponse = gfb_filter_zresponse(filter, z)
%
%   Input parameters:
%     filter  : A gfb_Filter struct as created by gfb_filter_new.
%     z       : A vector of z-plane frequencies where the frequency response should
%               be computed. $z = \exp(2i\cdot pi\cdot f/fs)$
%   Output parameters:
%     zresponse : The complex response of the filter at z.
%
%   Computes the frequency response of the gammatone filter at the frequency z.

% author   : tp
% date     : Jan & Nov 2006

warning('Warning: GFB_FILTER_ZRESPONSE will be removed in a future release. Use HOHMANN2002FREQZ instead. ');

zresponse = (1 - filter.coefficient ./ z) .^ -filter.gamma_order * ...
    filter.normalization_factor;
