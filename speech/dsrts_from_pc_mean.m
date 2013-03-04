function [dSRTs_out newSRT] = dSRTs_from_PC_mean(PC_input,SNRs,conditions)
%%
% dSRTs_from_PC_mean: calculates the SRT and change in SRT from the simlated percent correct 
% 
% Usage: [dSRTs_out newSRT] = dSRTs_from_PC_mean(PC_input,SNRs,conditions)
% PC_input   :  Matrix with the mean percent correct for each processing
%               condition and SNR. The first column of PC_input should always be the
%               reference with no processing
% conditions :  Vector with the processing conditions
% SNRs       :  Vector with the SNRs used.
% 
% Søren Jørgensen august 2010
%

 % ----------------  calculating dSRTs on mean psychofuncs.
for k = conditions
    SNRsLong = 0;
    yhat = 0;
    % ----------- Connecting the points with streight lines-----------
    for p = 1:length(SNRs)-1
        tmp =  polyfit(SNRs(p:p+1), PC_input(p:p+1,k)',1);
        tmpSNR =  SNRs(p):0.005:SNRs(p+1);
        yhat_tmp = polyval(tmp,tmpSNR);
        
        yhat = [yhat yhat_tmp];
        SNRsLong = [SNRsLong tmpSNR];
        
    end
    
    yhat2(:,k) = yhat(2:end);
    
    SNRsLong = SNRsLong(2:end);
    SRT = 50;
    hLine = linspace(SRT,SRT,length(yhat2(:,k)))';
    C = abs(yhat2(:,k)-hLine); %  Find minimum distance between horizontal line and estimated p_correct
    SNRIndex = find(C==min(C));
    if length(SNRIndex)>1
        SNRIndex = SNRIndex(1);
    end
    newSRT(k) =  SNRsLong(SNRIndex);
    
%     ----------- calculating change in SRT ---------------------
    dSRTs = abs(newSRT(1)-newSRT(k));
    if newSRT(k) > newSRT(1)
        dSRTs = dSRTs;
    else
        dSRTs = dSRTs*-1;
    end
    if min(C) > 1
        dSRTs = -4;
    end
    if  isempty(dSRTs)
        dSRTs = -4;
    end
    dSRTs_out(k) = dSRTs;
end
