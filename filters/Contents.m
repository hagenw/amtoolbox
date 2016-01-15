% AMT - Filter functions
%
%   The AMT team, 2012 - 2014.
%
%   General routines
%     UFILTERBANKZ     - Apply multiple filters
%     FILTERBANKZ      - Apply multiple filters with non-equidistant downsampling
%     FILTERBANK_INIT  - Create control structure for FILTERBANK_BLOCK
%     FILTERBANK_BLOCK - Filterbank block processing
%
%   Auditory filters
%     GAMMATONE        - Calculate Gammatone filter coefficients
%     CQDFT            - FFT-based filter bank with constant relative bandwidth
%
%   Hohmann (2002) Gammatone filterbank framework
%     GFB_ANALYZER_NEW            - Gammatone filterbank implementation, see demo_hohmann2012
%     GFB_ANALYZER_PROCESS        - All gfb functions are part of the hohmann2012 model
%     GFB_DELAY_NEW               - .
%     GFB_DELAY_PROCESS           - .
%     GFB_FILTER_NEW              - .
%     GFB_FILTER_PROCESS          - .
%     GFB_MIXER_NEW               - .
%     GFB_SYNTHESIZER_NEW         - .
%     HOHMANN2002CLEARSTATE       - Clears the state of the filterbank objects
%     HOHMANN2002FREQZ            - Calculates frequency response of a filterbank object
%
%   Averaging
%     WEIGHTEDAVERAGEFILTER       - Part of the takanen2013 model
%
%  For help, bug reports, suggestions etc. please send email to
%  amtoolbox-help@lists.sourceforge.net

% Copyright (C) 2009-2015 Piotr Majdak
% This file is part of AMToolbox version 0.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

