% This file defines global constants for the matlab gammatone filterbank
% implementation.
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2006

% filename : Gfb_set_constants

global GFB_L GFB_Q GFB_PREFERED_GAMMA_ORDER GFB_GAINCALC_ITERATIONS;

GFB_L = 24.7;  % see equation (17) in [Hohmann 2002]
GFB_Q = 9.265; % see equation (17) in [Hohmann 2002]

% We will use 4th order gammatone filters:
GFB_PREFERED_GAMMA_ORDER = 4;

% The gain factors are approximated in iterations. This is the default
% number of iterations:
GFB_GAINCALC_ITERATIONS  = 100;
