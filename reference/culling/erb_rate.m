function rate = erb_rate(f)
% f is in Hz, not kHz
f = f / 1000;
rate = 11.17 * log((f+0.312)/(f+14.675)) + 43;