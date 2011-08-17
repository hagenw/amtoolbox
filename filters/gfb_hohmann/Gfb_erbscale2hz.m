function Hz = Gfb_erbscale2hz(ERBscale)
% Hz = Gfb_erbscale2hz(ERBscale)
% 
% implements equation (17) of [Hohmann 2002]: computes a frequency
% in Hz from its value on the ERBscale
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_set_erbscale2hz.m


global GFB_L;
global GFB_Q;
Gfb_set_constants;

Hz = (exp(ERBscale / GFB_Q) - 1) * (GFB_L * GFB_Q);

