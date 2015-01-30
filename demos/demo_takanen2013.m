%DEMO_TAKANEN2013 Demo of the binaural model by Takanen, Santala and Pulkki
%
%   This script generates a figure showing the result of the binaural
%   auditory model by Takanen, Santala and Pulkki (2013) for sound source
%   distributions consisting of different number of sound sources simulated
%   with HRTFs to emit incoherent samples of pink noise. The resulting
%   activity map shows that the activation spreads as the width of the
%   distribution increases, which is in accordance with the results of the
%   psychoacoustical experiment by Santala and Pulkki (2011).
%
%   Optionally, pre-computed cochlear model outputs can be applied to 
%   significantly reduce the required computation time. The pre-computed 
%   cochlear model outputs can be obtained from the authors.
%
%   .. figure::
%
%      Output of the audiory model
%
%      The activity map.
%
%   See also: takanen2013
%
%   References: takanen2013 takanen2014 santala2011

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland

%% Starting of the script
% Use pre-computed cochlear model outputs, otherwise set preComp=0;
preComp = 1;

compType = 1;
printFigs = 0;
printMap = 1;

%if the user wishes to use pre-computed cochlea model outputs to reduce the
%required computation time
if preComp ==1
    filename='demo_cochleadata.mat';
    data=amtload('takanen2013',filename);
    output= takanen2013(data.tests.cochlea,data.tests.fs,compType,printFigs,printMap);
    title(data.tests.scenario);
    set(gca,'Ytick',data.tests.ytickPos);set(gca,'YtickLabel',data.tests.ytickLab(end:-1:1));
    ylabel(data.tests.ylab);
%otherwise, binaural input signals are used
else
    filename='demo_binsignals.mat';
    data=amtload('takanen2013',filename);
    output= takanen2013(data.tests.insig,data.tests.fs,compType,printFigs,printMap);
    title(data.tests.scenario);
    set(gca,'Ytick',data.tests.ytickPos);set(gca,'YtickLabel',data.tests.ytickLab(end:-1:1));
    ylabel(data.tests.ylab);
end