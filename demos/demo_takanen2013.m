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
%     This is the headline of the figure, one line at most
%
%     This is the description of the figure. It can be
%     several lines.
%
%   See also: takanen2013
%
%   References: takanen2013a santala2011

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland

%% Starting of the script
preComp = 1;

compType = 1;
printFigs = 0;
printMap = 1;

f=[amtbasepath,filesep,'demos',filesep,'demo_takanen2013'];

%if the user wishes to use pre-computed cochlea model outputs to reduce the
%required computation time
if preComp ==1
    s=load([f, 'cochleadata.mat']);
    tests=s.tests;
    output= takanen2013(tests.cochlea,tests.fs,compType,printFigs,printMap);
    title(tests.scenario);
    set(gca,'Ytick',tests.ytickPos);set(gca,'YtickLabel',tests.ytickLab(end:-1:1));
    ylabel(tests.ylab);
%otherwise, binaural input signals are used
else
   s=load([f,'binsignals.mat']);
   tests=s.tests;
   output= takanen2013(tests.insig,tests.fs,compType,printFigs,printMap);
   title(tests.scenario);
   set(gca,'Ytick',tests.ytickPos);set(gca,'YtickLabel',tests.ytickLab(end:-1:1));
   ylabel(tests.ylab);
end