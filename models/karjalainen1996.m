function [slow,fast]=karjalainen1996(insig,fs)
%KARJALAINEN1996  Non-linear adapation network
%   Usage: [slow,fast]=karjalainen1996(insig,fs)
%
%   `[slow,fast]=karjalainen1996(insig,fs)` computes the non-linear
%   adaptation network from Karjalainen et al. (1996) on the input signal
%   *insig* sampled with a sampling frequency of *fs* Hz.
% 
%   * XXX What are the assumptions on the input? The example below
%     generates NaN's
%
%   * XXX What is *slow* and *fast*? In which units are they defined?
%
%   * XXX Which level convention is used ?
%
%   * XXX Bibtex entry for the correct paper to cite.
%
%   Examples:
%   ---------
%
%   The following show the adapation to a simple test signal:::
%
%     [insig,fs] = greasy;
%     [slow,fast]=karjalainen1996(insig,fs);
%
%     subplot(1,2,1);
%     plot(slow);
%
%     subplot(1,2,2);
%     plot(fast);
%
%   This file (and the corresponding mex file) where originally part of
%   HUTear- Matlab toolbox for auditory modeling. The toolbox is available at 
%   `<http://www.acoustics.hut.fi/software/HUTear/>`_

% AUTHOR: Aki Härmä, Helsinki University of Technology, 

[slow,fast]=comp_karjalainen1996(insig,fs);