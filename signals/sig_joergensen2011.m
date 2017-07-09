function [ir,fs]=sig_joergensen2011(duration)
%sig_joergensen2011  Return a simulated room impulse response for Joergensen models
%   Usage: ir=sig_joergensen2011(duration);
%          [ir,fs]=sig_joergensen2011(duration);
%
%   `sig_joergensen2011(duration)` returns a simulated impulse
%   response of length *duration*, where the duration is measured in
%   seconds. Only four durations (T30) are possible to select: 0.4, 0.7,
%   1.3, and 2.3 seconds.
%
%   `[ir,fs]=sig_joergensen2011(duration)` also returns the
%   sampling frequency, *fs* = 44100 Hz.
%
%   The impulse responses were created using the ODEON room acoustic
%   software version 10 (Christensen, 2009). The simulated room was shaped
%   like an auditorium with maximal dimensions of 28x16x10m
%   (length-width-height).
%
%   The source and the receiver were horizontally with a fixed distance of
%   5m, and placed approximately the center of the room. All surfaces had
%   the same coefficient, which was adjusted individually frequency such
%   that the room had similar reverberation (T30) in the octave bands from
%   63 to 8000 Hz.
%
%   The corresponding clarity (C50), defined as the ratio of the energy of
%   first 50 ms of the impulse response to the energy of the part, was 0.60,
%   -2.9, -6.6, and -8.0 dB for the four different lengths, respectively.
%
%   See also: joergensen2011 joergensen2013

switch(duration)
 case {0.4,0.7,1.3,2.3}
   [ir,fs]=amt_load('signals',['sig_joergensen2011_',num2str(duration),'s.wav']);
 otherwise
   error('%s: Unsupported duration.',upper(mfilename));  
end;
