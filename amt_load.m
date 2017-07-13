function varargout=amt_load(model,data,variable)
%amt_load Load auxiliary data of a model
%   Usage: amt_load(MODEL, DATA);
%
%   `amt_load(model, data)` loads the auxiliary data from the file `data`. The data will loaded 
%   from the directory `model` located in the local `auxdata` directory given by
%   `amt_auxdatapath`. 
%
%   If the file is not in the local `auxdata` directory, it will be downloaded from
%   the web address given by `amt_auxdataurl`.
%
%   The following file types are supported:
%     `.wav`: output will be as that from `audioread`
%     `.mat`: output as that as from `load`
%     others: output is the absolute filename
%
%   `amt_load(model, data, variable)` loads just a particular `variable`
%   from the file. 
%
%   See also: amt_auxdatapath amt_auxdataurl
%

  
%   Author: Piotr Majdak, 2015

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