% script for getting the (DRNL) filter parameters (Lopez-Poveda, meddis 2001)
% Author: Morten Løve Jepsen, 2.nov 2005, rev. 15 feb 2006, 19 feb 2007
%
% usage:  [linDRNLparOut,nlinDRNLparOut] = getDRNLparam(CF);

 function [linDRNLparOut,NlinDRNLparOut] = getDRNLparam(CF);

 %DRNL for normal hearing, Morten 2007

 % init structure
 linDRNLstruct  = struct('parname',{},'vals',{}); % initialize linear paramater vector
 NlinDRNLstruct = struct('parname',{},'vals',{}); % initialize nonlinear paramater vector
 linDRNLparOut  = struct('parname',{},'vals',{});  NlinDRNLparOut = struct('parname',{},'vals',{});
 
 linDRNLstruct(1).parname = 'CF_lin';  linDRNLstruct(2).parname = 'nGTfilt_lin';
 linDRNLstruct(3).parname = 'BW_lin';  linDRNLstruct(4).parname = 'g'; 
 linDRNLstruct(5).parname = 'LP_lin_cutoff';  linDRNLstruct(6).parname = 'nLPfilt_lin'; 
 linDRNLparOut=linDRNLstruct;
 
 NlinDRNLstruct(1).parname = 'CF_nlin';  NlinDRNLstruct(2).parname = 'nGTfilt_nlin'; 
 NlinDRNLstruct(3).parname = 'BW_nlin';  NlinDRNLstruct(4).parname = 'a'; 
 NlinDRNLstruct(5).parname = 'b';  NlinDRNLstruct(6).parname = 'c'; 
 NlinDRNLstruct(7).parname = 'LP_nlin_cutoff';  NlinDRNLstruct(8).parname = 'nLPfilt_nlin'; 
 NlinDRNLparOut=NlinDRNLstruct;

linDRNLstruct(1).vals = 10^(-0.06762+1.01679*log10(CF)); % Hz, CF_lin, 
linDRNLstruct(2).vals = 2; % number of cascaded gammatone filters 
linDRNLstruct(5).vals = 10^(-0.06762+1.01*log10(CF)); % Hz, LP_lin cutoff
    linDRNLstruct(6).vals = 4; % no. of cascaded LP filters, orig 4
    NlinDRNLstruct(1).vals = 10^(-0.05252+1.01650*log10(CF)); % Hz, CF_nlin
    NlinDRNLstruct(2).vals = 2; % number of cascaded gammatone filters 
    % In some code, this parameter has been .7, but in the jepsen2008
    % this parameter is .77
    NlinDRNLstruct(3).vals = 10^(-0.03193+.77*log10(CF)); % Hz, BW_nlin
    NlinDRNLstruct(7).vals = 10^(-0.05252+1.01*log10(CF)); % LP_nlincutoff
    NlinDRNLstruct(8).vals = 1; % no. of cascaded LP filters in nlin path, use "1"

    linDRNLstruct(4).vals = 10^(4.20405 -.47909*log10(CF)); %g
    linDRNLstruct(3).vals = 10^(.03728+.75*log10(CF)); % Hz, BW_lin. 
    if CF<=1000
        NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF)); % a, the 1500 assumption is no good for compressionat low freq filters
        NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(CF)); % b [(m/s)^(1-c)]
    else
        NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500)); % a, the 1500 assumption is no good for compressionat low freq filters
        NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(1500)); % b [(m/s)^(1-c)]
    end
    NlinDRNLstruct(6).vals = 10^(-.60206); % c, compression coeff

    for k=1:6
        linDRNLparOut(k).vals = linDRNLstruct(k).vals;
    end
    for k=1:8
        NlinDRNLparOut(k).vals = NlinDRNLstruct(k).vals;
    end
 
%OLDFORMAT
