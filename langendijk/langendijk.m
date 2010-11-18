function [ p ] = langendijk( ir1,ir2,bw,do,cp,s,bal,fstart,fend,stim )
% LANGENDIJK Localization model according to Langendijk et al. (2002)
% Usage:    [ p ] = langendijk( ir1,ir2 )
%           [ p ] = langendijk( ir1,ir2,stim )
%           [ p ] = langendijk( ir1,ir2,bw,do,cp,s,bal,fstart,fend )
%           [ p ] = langendijk( ir1,ir2,bw,do,cp,s,bal,fstart,fend,stim )
% Input arguments:
%     ir1:      modified impulse responses of DTFs for all positions and
%               both ears (sorted)
%     ir2:      stored impulse responses of DTF-templates for all positions 
%               and both ears (sorted)
%     bw:       bandwidth of filter bands as partial of an octave; default: 6
%     do:       differential order; default: 0
%     cp:       comparison process; 'std' (default) or 'xcorr'
%     s:        standard deviation of transforming Gaussian function; default: 2
%     bal:      balance of left to right channel (not included in 
%               langendijk's original comparison process); default: 1
%     fstart:   start frequency of filter bank; min: 0,5kHz; default: 2kHz
%     fend:     end frequency of filter bank; default: 16kHz
%     stim:     applied stimulus for localization test (optional)
% Output argument:
%     p:        predicted probability density functions for response angles 
%               with respect to target positions
% Description:
% p = langendijk( ir1,ir2,... ) results to a two dimensional matrix p. 
% The first dimension represents all possible response positions in 
% increasing order and the second dimension all possible target 
% respectively source positions. Consequently each column describes the 
% probability density function of the response distribution for one special 
% target position.
% If you want to plot this matrix use plotlangendijk().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-07-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default settings according to Langendijk et al. (2002)
if ~exist('bw','var')
    bw=6;
end
if ~exist('do','var')
    do=0;
end
if ~exist('cp','var')
    cp='std';
end
if ~exist('s','var')
    s=2;
end
if ~exist('bal','var')
    bal=1;
end
if ~exist('fstart','var')
    fstart=2000;
end
if ~exist('fend','var')
    fend=16000;
end

% addition for applied stimulus
if length(bw)>1
    stim=bw;
end
if exist('stim','var')
    ir1=halfconv( ir1,stim );
end

% model calculations
p=zeros(size(ir2,2),size(ir1,2)); % initialisation
% filter bank
for ind2=1:size(ir1,2) % response pdf for every target position
    x=averagingfb(ir1(:,ind2,:),bw,fstart,fend);
    y=zeros(length(x),size(ir2,2),size(ir2,3)); % initialisation
    for ind=1:size(y,2) % response pdf for one target position
        y(:,ind,:)=averagingfb(ir2(:,ind,:),bw,fstart,fend);
    end
% comparison process for one target position
    p(:,ind2)=langecomp( x,y,s,do,cp,bal );
end
end