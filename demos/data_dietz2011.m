function [stim,fs]  = data_dietz2011(varargin)
%DATA_DIETZ2011 Stimuli from Dietz et al. (2011)
%   Usage: [stim,fs] = data_dietz2011;
%
%   `[stim,fs]=data_dietz2011(flag)` returns the speech stimuli and sampling
%   frequency from Dietz et al. 2011. The possible values of *flag* are:
%   'five_speakers', 'one_of_three', 'one_speaker_reverb', 'three_of_three',
%   'bNoise', 'two_of_three', or 'two_speakers'.
%
%   See also: dietz2011, exp_dietz2011
%
%   References: dietz2011auditory

definput.flags.wavfile={'missingflag','five_speakers','one_of_three',...
                    'one_speaker_reverb','three_of_three',...
                    'two_of_three','two_speakers','bNoise'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.wavfile{2:end-2}),...
               sprintf('%s or %s',definput.flags.wavfile{end-1},definput.flags.wavfile{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;
  
s=[amtbasepath,'humandata',filesep,'dietz2011_',flags.wavfile,'.wav'];

[stim,fs] = wavread(s);

  