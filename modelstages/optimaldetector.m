function mue = optimaldetector(ir_stim,template)
%OPTIMALDETECTOR  Generic optimal detector for the CASP and Breebaart models
%
%   This is a correlation-based optimal detector for a signal known exactly
%   see Green & Swets (1966) for more details.

corrmue = ir_stim.*template;
optfactor = sqrt(numel(corrmue));

% Take mean over all dimensions of internal representation and correct for
% optimalityfactor.
mue = mean(corrmue(:))*optfactor;

