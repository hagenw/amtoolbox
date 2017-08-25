function [spikeRatePerNeuron,spikeRatePerBin] = ...
                        kelvasa2015anbinning(APvec,sigLengthSec,...
                                                    varargin)
%KELVASA2015ANBINNING  AN and time binning from Kelvasa and Dietz 2015 binaural model
%   Usage: [spikeRatePerNeuron,spikeRatePerBin] = ...
%                        kelvasa2015anbinning(APvec,sigLengthSec);
%
%   Input parameters:
%       APvec         : N x 2 matrix of AN spikes with Nx1 holding indices 
%                       of the spiking neuron and Nx2 holding corresponding
%                       spike time in seconds. 
%
%       sigLengthSec  : length of input signal in seconds
%
%   Output parameters:
%       spikeRatePerNeuron: N x M matrix of AN spike rates in spikes/second
%                           with N being the number of user defined AN 
%                           fibers and M being the number of time windows.
%
%       spikeRatePerBin   : N x M matrix of AN spike rates in spikes/second
%                           with N being the number of user defined AN fibe
%                           bands and M being the number of time windows. 
%
%   KELVASA2015anbinning(APvec,sigLengthSec,varargin) bins auditory nerve 
%   spike times over a given population of AN fibers into user defined AN 
%   frequency bands and time bins as detailed in (Kelvasa & Dietz (2015))
%
%   References: kelvasa2015
         
%   Authors: 
%            Daryl Kelvasa (daryl.kelvasa@uni-oldenburg.de) 2016
%            Mathias Dietz (mdietz@uwo.ca) 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input paramters
if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

%Retrieve and compute model paramters
definput.import={'kelvasa2015'};
[~,kv]  = ltfatarghelper({},definput,varargin);
                    
%% Main code

%initialize variables
numWin = ceil(sigLengthSec/kv.timeWindowSec);
winEdges = linspace(0,sigLengthSec,numWin+1);
spikeRatePerNeuron = zeros(kv.N_nervecells,numWin);
spikeRatePerBin = zeros(kv.numBin,numWin);

if ~isempty(APvec)
    
    [~,ind] = histc(APvec(:,2),winEdges);
    ind(ind==numWin+1) = numWin;
    
    for win = 1 : numWin
   
        APwin = APvec(ind == win,:);
        [H, ~] = histc(APwin(:,1),1:kv.N_nervecells);
        clear APwin
        spkRate =  H./kv.timeWindowSec;
        spikeRatePerNeuron(:,win) = spkRate;
        spikeRatePerBin(:,win) = mean(spkRate(kv.binPosInd),2);
    
    end
end

end
