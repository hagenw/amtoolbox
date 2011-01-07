function decision=breebart2001(intervals,fs,tau,alpha)
%BREEBART2001  The Breebart2001 model
%
%  Input arguments:
%    Intervals  - Maxtrix of intervals. Dimensions: time, channel,
%                 interval no.
%
%
%  
  
% Compute the preprocessing of all the intervals
ir_all = breebart2001preproc(intervals,fs);

[siglen,nfreqchannels,naudiochannels,nifc] = size(ir_all);

ei_map = zeros(nifc, nfreqchannels, siglen);
for k=1:nifc
  for g=1:nfreqchannels
    ei_map(k,g,:) = eicell(squeuze(ir(:,:,g,k)),fs,tau,alpha);
  end
end


