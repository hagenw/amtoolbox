function synth = hohmann2002synth(fb, desired_delay)
%HOHMANN2002SYNTH  Create new synthesizer object within HOHMANN2002 filterbank framework
%   Usage: synth = hohmann2002synth(fb, desired_delay)
%
%   Input parameters:
%     fb            : The filterbank structure as returned by |hohmann2002|.
%     desired_delay : the desired group delay of the total analysis-synthesis
%                     system (in seconds). Minimum delay is 1 sample. 
%                     Greater delays result in better output signal quality. 
%
%   Output parameters:
%     synth : the constructed synthesizer structure
%
%   `hohmann2002synth` creates a new synthesizer object for the
%   reconstruction of signals analyzed by `fb` within the HOHMANN2002
%   filterbank framework.

% author: Universitaet Oldenburg, tp (Jan 2002, Jan, Sep 2003, Nov 2006, Jan 2007)
% Adapted to AMT (PM, Jan 2016) from function gfb_synthesizer_new

synth.type = 'gfb_Synthesizer';
desired_delay_in_samples = round(desired_delay * fb.fs);
if (desired_delay_in_samples < 1)
    error('delay must be at least 1 sample');
end

synth.delay = hohmann2002delay(fb, desired_delay_in_samples);
synth.mixer = hohmann2002mixer(fb, synth.delay);
