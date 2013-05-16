function [ si ] = langendijk2002comp( target,template,varargin )
%LANGENDIJK2002COMP Comparison process according to Langendijk et al. (2002)
%   Usage:        [ si ] = langendijk2002comp( target,template )
%                 [ si ] = langendijk2002comp( target,template,s,do,flags )
%
%   Input parameters:
%     target   :  modified DFT for one target position and both ears
%     template :  internal DTF-templates for all positions and both ears
%
%   Output parameters:
%     si       :  monaural similarity indices (one for each template entry)
%
%   `langendijk2002comp` compares the spectral features of an internal spectral
%   representation of a target sound with an internal template. As spectral
%   similarity criterion either the crosscorrelation coefficient ('xcorr') 
%   or a mapped version of the standard deviation of inter-spectral 
%   differences ('std') can be used.
%
%   The function accepts the following optional parameters:
%
%     's', s     Standard deviation of transforming Gaussian function 
%                (only for comparison process: std).
%                The default value is 2.
%
%     'do',do    Differential order. The default value is 0.
%
%     'std'      Use the 'std' comparison process. This is the default.
%
%     'xcorr'    Use the 'xcorr' comparison process. 
%
%   See also: langendijk2002
  
  
% AUTHOR: Robert Baumgartner, OEAW Acoustical Research Institute


definput.import={'langendijk2002comp'};
[flags,kv]=ltfatarghelper({'s','do'},definput,varargin);


%% Optional differentiation

if kv.do ~= 0
	target = diff(target,kv.do);
    template = diff(template,kv.do);
end


%% Comparison process

if flags.do_std
    
    idiff = repmat(target,[1,size(template,2),1]) - template;
    sigma = squeeze(std(idiff,1));
    si = normpdf(sigma,0,kv.s);
    
elseif flags.do_xcorr
    
    r = zeros(size(template,2),size(template,3));
    for ch = 1:size(template,3)
        for ii = 1:size(template,2)
            tmp = corrcoef(target(:,ch),template(:,ii,ch));
            r(ii,ch) = tmp(2);
        end
    end
    si = (r+1)/2;
    
end

end



% if flags.do_std
%   
%   ptemp=zeros(size(in2,2),size(in2,3)); % initialisation
%   for ch=1:size(in2,3)
%     for ind=1:size(in2,2)
%       if kv.do==0
%         z=in1(:,ch)-in2(:,ind,ch);
%       else
%         z=diff(in1(:,ch),kv.do)-diff(in2(:,ind,ch),kv.do);
%       end
%       sigma=std(z,1);
%       ptemp(ind,ch)=normpdf(sigma,0,kv.s);
%     end
%   end
%   if size(in2,3)==1
%     p=ptemp;
%   else
%     p=(kv.bal*ptemp(:,1)+1/kv.bal*ptemp(:,2))/2; % balance
%   end
%   p=1/sum(p).*p;
%   % p=p/max(p); % normalisation
%     
% end;
%   
% if flags.do_xcorr
% 
%   ptemp=zeros(size(in2,2),size(in2,3)); % initialisation
%   for ch=1:size(in2,3)
%     for ind=1:size(in2,2)
%       if kv.do==0
%         z=corrcoef(in1(:,ch),in2(:,ind,ch));
%       else
%         z=corrcoef(diff(in1(:,ch),kv.do),diff(in2(:,ind,ch),kv.do));
%       end
%       z=z(2);
%       ptemp(ind,ch)=(z+1)/2;
%     end
%   end
%   if size(in2,3)==1
%     p=ptemp;
%   else
%     p=(kv.bal*ptemp(:,1)+1/kv.bal*ptemp(:,2))/2; % balance
%   end
%   
% end
% 
% out=p;


