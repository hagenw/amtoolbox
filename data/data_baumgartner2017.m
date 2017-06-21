function data = data_baumgartner2017(varargin)
%DATA_BAUMGARTNER2017  Data from Baumgartner et al. (2017)
%   Usage: data = data_baumgartner2017(flag)
%
%   `data_baumgartner2017(flag)` returns data from Baumgartner et al. (2017)
%   describing a model for sound externalization.
%
%   The input flag determines the microphone casing and may be
%
%     'ITE'       to obtain a set of 23 standard in-the-ear HRTFs (default)
%     'BTE'       to obtain an exemplary behind-the-ear HRTF
%
%   The fields in the *data* output contains the following information
%
%     .id         listener ID
%     .Obj        HRTF data in SOFA Format
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2017

% AUTHOR : Robert Baumgartner

% TODO: explain Data in description;

%% % Parse input options
definput.flags.type = {'ITE','BTE'};
flags = ltfatarghelper({},definput,varargin);
   
dataPath = fullfile(SOFAdbPath,'baumgartner2017');
if flags.do_ITE
  fn = dir(fullfile(dataPath,'hrtf b_nh*.sofa'));
elseif flags.do_BTE
  fn = dir(fullfile(dataPath,'hrtf b_bte_nh*.sofa'));
end
if isempty(fn)
  error('RB: No files found. Check directories and AMT setup!')
end
data = struct('id',{},'Obj',{});
for ii = 1:length(fn)
  data(ii).id = fn(ii).name(end-6:end-5);
  data(ii).Obj = SOFAload(fullfile(dataPath,fn(ii).name));
end

end