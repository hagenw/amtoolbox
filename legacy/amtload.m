function varargout=amtload(model,data,variable)

warning('Warning: AMTLOAD will be removed in a future release. Use AMT_LOAD instead. ');  

% based on amt_load from 9.7.2017

localfn=fullfile(amt_auxdatapath,model,data);
  % file not found? create directories, and download!
if ~exist(localfn,'file')
    % create dir if not existing
  if ~exist(fullfile(amt_auxdatapath,model),'dir'), 
    [success,msg]=mkdir(fullfile(amt_auxdatapath,model));
    if success~=1
      error(msg);
    end
  end
    % download
  amt_disp(['Model: ' model '. Downloading auxiliary data: ' data],'progress');
  webfn=[amt_auxdataurl '/' model '/' data];
  webfn(strfind(webfn,'\'))='/';
  webfn=regexprep(webfn,' ','%20');        
  [~,stat]=urlwrite(webfn,localfn);
  if ~stat
    error(['Unable to download file: ' webfn]);
  end
end
  % load the content
[~,~,ext] = fileparts(localfn);
switch lower(ext)
  case '.wav'
    [y,fs]=audioread(localfn);
    varargout{1}=y;
    varargout{2}=fs;
  case '.mat'
    if exist('variable','var'),
        varargout{1}=load(localfn,variable);
    else
        varargout{1}=load(localfn);
    end
  otherwise
    varargout{1}=localfn;
end