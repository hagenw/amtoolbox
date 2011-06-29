function [medir]=gr2ir(med,cond,fs)
% GR2IR converts given gain responses MED (in dB) to impulse responses
% MEDIR; furthermore several conditions according to langendijk et al.
% (2002) can be defined
% Usage:            [medir]=gr2ir(med,cond,fs)
% Input arguments:
%       med:        gain responses (in dB)
%       cond:       condition, 
%                   possibilities:  baseline    'b'
%                                   2 octaves   '2o'
%                                   1 oct (low) '1ol'
%                                   1 oct (mid) '1om'
%                                   1 oct (high)'1oh'
%       fs:      	sampling frequency
% Output arguments:
%       medir:   	impulse responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default settings
if ~exist('cond','var')
    cond='b';
end
if ~exist('fs','var')
    fs=48000;
end

imp=zeros(240,1);imp(1)=1;
n=256;
len=240;
frq=logspace(log10(2000),log10(16000),size(med,1));

% frequency indices
f1=find(frq>=4000,1);
f2=find(frq>=5700,1);
f3=find(frq>=8000,1);
f4=find(frq>=11300,1);

switch cond
    case 'b' % baseline
        medir=zeros(n,size(med,2));
        for ii=1:size(med,2)
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '2o' % 2 oct (4-16kHz)
        medir=zeros(n,size(med,2));
        med2o=med;
        for ii=1:size(med,2)
            med2o(f1:end,ii)=mean(med(f1:end,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med2o(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1ol' % 1 oct (low:4-8kHz)
        medir=zeros(n,size(med,2));
        med1ol=med;
        for ii=1:size(med,2)
            med1ol(f1:f3,ii)=mean(med(f1:f3,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1ol(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1om' % 1 oct (middle:5.7-11.3kHz)
        medir=zeros(n,size(med,2));
        med1om=med;
        for ii=1:size(med,2)
            med1om(f2:f4,ii)=mean(med(f2:f4,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1om(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end

    case '1oh' % 1 oct (high:8-16kHz)
        medir=zeros(n,size(med,2));
        med1oh=med;
        for ii=1:size(med,2)
            med1oh(f3:end,ii)=mean(med(f3:end,ii));
            bp = firls(len,[0 0.08 frq/(fs/2) 0.68 1],[0;0; 10.^(med1oh(:,ii)/20); 0;0]);
            medir(1:length(imp),ii) = filter(bp,1,imp);
        end
end
end
