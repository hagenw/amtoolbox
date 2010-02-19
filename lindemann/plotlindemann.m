function plotlindemann(crosscorr,fc)

% ------ Checking of input  parameters -----------------------------------

error(nargchk(1,2,nargin));

if ~isnumeric(crosscorr)
    error('%s: crosscorr has to be numeric!',upper(mfilename));
end

% ------ Computation -----------------------------------------------------
    
% Calculate tau (delay line time)
tau = linspace(-1,1,size(crosscorr,2));
% FIXME: calculate t
t = 0:size(crosscorr,1)-1;
% Calculate mean binaural activation pattern
if nargin==1 || ~isscalar(fc)
    binpattern = mean(crosscorr(1:end,:,:),3);
else
    binpattern = crosscorr(1:end,:,fc);
end

% plot result
figure;
mesh(tau,t,binpattern);