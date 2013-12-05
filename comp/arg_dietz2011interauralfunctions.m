function definput=arg_dietz2011interauralfunctions(definput)

  definput.keyvals.signal_level_dB_SPL = 70;
  definput.keyvals.tau_cycles  = 5; % see Fig. 3c in Dietz (2011)
  definput.keyvals.compression_power = 0.4;

  % ask for simulating temporal resolution of binaural processor
  % this returns the *_lp values
  definput.flags.lowpass = {'lowpass','nolowpass'};
