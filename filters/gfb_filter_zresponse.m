function zresponse = gfb_filter_zresponse(filter, z)
% zresponse = gfb_filter_zresponse(filter, z)
%
% Computes the frequency response of the gammatone filter at the frequency z.
% 
%  Input parameters:
% filter  A gfb_Filter struct as created by gfb_filter_new.
% z       A vector of z-plane frequencies where the frequency response should
%         be computed. z = exp(2i*pi*f[Hz]/fs[Hz])
% zresponse The complex response of the filter at z.
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan & Nov 2006

zresponse = (1 - filter.coefficient ./ z) .^ -filter.gamma_order * ...
    filter.normalization_factor;

%OLDFORMAT
