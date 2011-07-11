function [template,ir_reference]=casptemplate(target,reference,modelname,modelpars)
%CASPTEMPLATE  Generate a template for the optimal detector
%
%  CASPTEMPLATE(target,reference,modelname,modelpars) generates the template
%  needed for the optimal detector. CASPTEMPLATE will run the model
%  specifief by modelname on the signals stored in target and reference
%  and generate the template from this.
%
%  If target or reference is a matrix, each column will be considered a
%  signal, and averaging will be done. This is usefull for stochastic
%  signals.

if nargin<4
  modelpars={};
end;

ntargets    = size(target,2);
nreferences = size(reference,2);

%% ----- Compute average internal representation of the targets
ir_target=feval(modelname,target(:,1),modelpars{:});

for ii=2:ntargets
  ir_target = ir_target + feval(modelname,target(:,ii),modelpars{:});
end;

ir_target=ir_target/ntargets;

%% ----- Compute average internal representation of the references
ir_reference=feval(modelname,reference(:,1),modelpars{:});

for ii=2:nreferences
  ir_reference = ir_reference + feval(modelname,reference(:,ii),modelpars{:});
end;

ir_reference=ir_reference/nreferences;

% Compute the template as the difference between the average
% representation of the targets and references.
template = ir_target - ir_reference;

%% ----- Normalize to compensate for the increase in level ----

% Normalize across all dimenstions of the internal representation.
template=template/rms(tempate(:));