function [ X ] = discreteinvrnd(p,m,n)
% DISCRETEINVRND implements an inversion method for a discrete distribution
% with probability mass vector p and dimensions m,n
% Usage:    [ X ] = discreteinvrnd(p,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW
% latest update: 2010-07-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(m,n); 
for i = 1:m*n
    c = cumsum(p);
    u = max(c)*rand;
    X(i) = find(u < c ,1);
end
end