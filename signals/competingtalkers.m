function [s,fs]=competingtalkers(varargin)
%COMPETINGTALKERS  Load one of several test signals
%   Usage:  s=competingtalkers(signame);
%           [s,fs]=competingtalkers(signame);
%
%   `competingtalkers(signame)` loads one of several test signals consisting
%   of competing talkers. All the talkers are taken from the TIMIT speech
%   corpus:
%   `<http://www.ldc.upenn.edu/Catalog/CatalogEntry.jsp?catalogId=LDC93S1>`_.
%
%   The signals have 2 channels and are all recorded with a sampling rate of
%   16 kHz.
%
%   `[sig,fs]=competingtalkers(signame)` additionally returns the sampling
%   frequency *fs*.
%
%   The value of `signame` may be one of:
%
%     'one_of_three'    XXX Description missing
%
%     'two_of_three'    XXX Description missing
%
%     'three_of_three'  XXX Description missing
%
%     'one_speaker_reverb'
%                       XXX Description missing
%
%     'two_speakers'    XXX Description missing
%
%     'five_speakers'   XXX Description missing
%
%     'bnoise'          Speech shaped noise
%
%   Examples:
%   ---------
%
%   The following plot shows an estimate of the power spectral density of
%   the first channels of the speech shaped noise:::
%
%      s=competingtalkers('bnoise');
%      pwelch(s(:,1),hamming(150));
%
%   See also: exp_dietz2011

%   AUTHOR : Peter L. Søndergaard


definput.flags.sigtype={'missingflag','one_of_three','two_of_three',...
                    'three_of_three','one_speaker_reverb',...
                    'two_speakers','five_speakers','bnoise'};

[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.sigtype{2:end-2}),...
             sprintf('%s or %s',definput.flags.sigtype{end-1},definput.flags.sigtype{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


% f=mfilename('fullpath');
% 
% fname=[f,'_',flags.sigtype,'.wav'];
% s = wavread(fname);
% fs = 16000;
[s,fs]=amtload('competingtalkers',[flags.sigtype '.wav']);