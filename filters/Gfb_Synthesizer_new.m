function synthesizer = Gfb_Synthesizer_new(analyzer, desired_delay_in_seconds)
% synthesizer = Gfb_Synthesizer_new(analyzer, desired_delay_in_seconds)
%
% Gfb_Synthesizer_new creates a new synthesizer object that fits to the
% given analyzer.
%
% PARAMETERS:
% analyzer                  an analyzer struct as returned by Gfb_Analyzer_new
% desired_delay_in_seconds  the desired group delay of the analysis-synthesis
%                           system in seconds.  Greater delays result in better
%                           output signal quality.  Minimum delay is
%                           (1 / analyzer.fs)
% synthesizer               the constructed Gfb_Synthesizer structure
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Nov 2003, Nov 2006, Jan 2007

% filename : Gfb_Synthesizer_new.m


synthesizer.type         = 'Gfb_Synthesizer';
desired_delay_in_samples = round(desired_delay_in_seconds * ...
			         analyzer.fs);
if (desired_delay_in_samples < 1)
    error('delay must be at least 1/analyzer.fs');
end

synthesizer.delay = Gfb_Delay_new(analyzer, desired_delay_in_samples);
synthesizer.mixer = Gfb_Mixer_new(analyzer, synthesizer.delay);

%OLDFORMAT
