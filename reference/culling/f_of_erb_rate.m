function f = f_of_erb_rate(erb)
f = (0.312 - (exp((erb - 43)/11.17)) * 14.675) / (exp((erb - 43)/11.17) - 1);
f = f * 1000; %convert f to Hz not kHz