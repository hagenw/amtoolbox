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
%   Gammatone filterbank framework from Hohmann (2002)
%     HOHMANN2002DELAY            - Create a delay object
%     HOHMANN2002MIXER            - Create a mixer object
%     HOHMANN2002SYNTH            - Create a synthesis object 
%     HOHMANN2002FILT             - Create a single Gammatone filter object
%     HOHMANN2002PROCESS          - Process the input signals by the corresponding filterbank object
%     HOHMANN2002CLEARSTATE       - Clears the state of the filterbank object
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

