function ERBscale = Gfb_hz2erbscale(Hz)
% ERBscale = Gfb_hz2erbscale(Hz)
% 
% implements equation (16) of [Hohmann 2002]:  computes an ERBscale
% value from a frequency in Hz
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_hz2erbscale.m


global GFB_L;
global GFB_Q;
Gfb_set_constants;

ERBscale = GFB_Q * log(1 + Hz / (GFB_L * GFB_Q));

