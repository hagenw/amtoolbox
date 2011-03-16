function [t,table]=absolutethreshold(freq,varargin)
%ABSOLUTETHRESHOLD  Absolute threshold at specified frequencies
%   Usage:  t=absolutethreshold(freq);
%           t=absolutethreshold(freq,...);
%
%   ABSOLUTETHRESHOLD(freq) will return the absolute threshold of hearing in
%   dB SPL at the frequencies given in freq. The output will have the same
%   shape as the input.
%
%   [t,x,y]=ABSOLUTETHRESHOLD(...) additionally returns the frequencies x
%   and levels y specifying the choosen standard.
%
%   ABSOLUTETHRESHOLD takes the following optional parameters:
%
%-     'iso226_2003' - Free field binaural as given by the ISO 226:2003 standard.
%
%-     'map'         - The ISO 226:2003 standard coverted to minimal
%                      audible pressure using the method from Bentler &
%                      Pavlovic 1989.
%
%-     'er3a'        - Using ER-3A insert earphones. This is described in
%                      the ISO 389-2:1994(E) standard
%
%-     'er2a'        - Using ER-2A insert earphones. This is described in
%                      Han & Poulsen (1998)
%
%-     'hda200'      - Using HDA200 circumaural earphone. This is
%                      decribed in the ISO 389-8:2004 standard.
%
%   Absolute thresholds for the ER2A and HDA-200 are provided up to 16
%   kHz by the ISO 389-5:2006 standard.
%
%   The default is to use the 'iso226_2003' setting.
%
%   Demos: demo_absolutethreshold
%
%R  iso226-2003 iso389-2-1994 iso389-5-2006 iso389-8-2004 bentler1989transfer han1998equivalent
  
% AUTHOR : Peter Søndergaard based on data collected by Claus Elberling
  
if nargin<1
  error('Too few input parameters.');
end;

definput.import = {'absolutethreshold'};
definput.importdefaults = {'iso226_2003'};

[flags,kv,table]  = ltfatarghelper({'table'},definput,varargin);

t=spline(table(:,1),table(:,2),freq);
