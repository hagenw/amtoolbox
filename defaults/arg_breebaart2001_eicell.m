function definput=arg_breebaart2001_eicell(definput)

  definput.keyvals.tc    = 30e-3;   % Temporal smoothing constant
  definput.keyvals.rc_a  = 0.1;     % Range compression parameter 'a' 
  definput.keyvals.rc_b  = .00002;  % Range compression parameter 'b'
  definput.keyvals.ptau  = 2.2e-3;  % time constant for p(tau) function
