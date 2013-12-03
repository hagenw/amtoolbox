function definput=arg_dietz2011modulationfilter(definput)
 
  % Parameters for filtering the haircell output
  definput.keyvals.filter_order = 2;            % used for both env and fine
  definput.keyvals.filter_attenuation_db = 10;  % used for both env and fine
  % Finestructure filter
  definput.keyvals.fine_filter_finesse = 3;
  % Envelope/modulation filter
  definput.keyvals.mod_center_frequency_hz = 135;
  definput.keyvals.mod_filter_finesse = 8; % => bandwidth: 16.9 Hz
  % ILD filter
  definput.keyvals.level_filter_order = 2;
  definput.keyvals.level_filter_cutoff_hz = 30;

