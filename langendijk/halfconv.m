function [ out ] = halfconv( ir1,stim )
% HALFCONV calculates the fast convolution with fft but without ifft
% Usage:        [ out ] = halfconv( ir1,stim )
% Input arguments:
%     ir1:      (modified) impulse responses of DFTs for all positions and
%               both ears
%     stim:     stimulus (time domain)
% Output argument:
%     out:      convolution of stimulus with DTFs in frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robert Baumgartner, OEAW
% latest update: 2010-07-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfft = 2^nextpow2(max([size(ir1,1) length(stim)]));
stimf=fft(stim(:),nfft);
ir1f = fft(ir1,nfft,1);
temp=zeros(nfft,size(ir1,2),size(ir1,3));
for ch=1:size(ir1,3)
    for ind=1:size(ir1,2)
        temp(:,ind,ch) = ir1f(:,ind,ch).* stimf;
    end
end
out=temp;

end

