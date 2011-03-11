function [t,x,y]=absolutethreshold(freq,varargin)
%ABSOLUTETHRESHOLD  Absolute threshold at specified frequencies
%   Usage:  t=absolutethreshold(freq);
%           t=absolutethreshold(freq,...);
%
%   ABSOLUTETHRESHOLD(freq) will return the absolute threshold of hearing in
%   dB SPL at the frequencies given in freq. The output will have the same
%   shape as the input.
%
%   ABSOLUTETHRESHOLD takes the following optional parameters:
%
%-     'iso226_2003' - Free field binaural as given by the ISO 226:2003 standard.
%
%-     'map'         - The ISO 226:2003 standard coverted to monaural ear drum
%                      pressure using the method from Bentler &
%                      Pavlovic 1989.
%
%-     'er3a'        - Using ER-3A insert earphones. This is described in
%                      the ISO 389-2:1994(E) standard
%
%-     'er2a'        - Using ER-2A insert earphones. This is described in
%                      Han & Poulsen (1998)
%
%   The default is to use the 'iso226_2003' setting.
%
%R  bentler1989transfer han1998equivalent
  
% AUTHOR : Peter Søndergaard based on data collected by Claus Elberling
  
if nargin<1
  error('Too few input parameters.');
end;

definput.import = {'absolutethreshold'};
definput.importdefaults = {'iso226_2003'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

x=kv.deffreq;
y=kv.level;

t=spline(kv.deffreq,x,y);
  
  