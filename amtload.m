function varargout=amtload(model,data,variable)
%AMTLOAD Load auxiliary data of a model
%   Usage: amtload(MODEL, DATA);
%
%   `amtload` loads the auxiliary data from the file `data`. The data will loaded 
%   from the directory `model` located in the auxdata directory given by
%   `amtauxdatapath`. 
%
%   If the file is not in the auxdata directory, it will be downloaded from
%   the web address given by `amtauxdataurl`
%
%   The following file types are supported:
%     `.wav`: output will be as that from `audioread`
%     `.mat`: output as that as from `load`
%     others: output is the absolute filename
%
%
%   See also: amtauxdatapath amtauxdataurl
%

  
%   Author: Piotr Majdak, 2015

localfn=fullfile(amtauxdatapath,model,data);
  % file not found? create directories, and download!
if ~exist(localfn,'file')
    % create dir if not existing
  if ~exist(fullfile(amtauxdatapath,model),'dir'), 
    [success,msg]=mkdir(fullfile(amtauxdatapath,model));
    if success~=1
      error(msg);
    end
  end
    % download
  amtdisp(['Model: ' model '. Downloading auxiliary data: ' data],'progress');
  webfn=[amtauxdataurl '/' model '/' data];
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