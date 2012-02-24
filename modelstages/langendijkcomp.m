function [ out ] = langendijkcomp( in1,in2,varargin )
%LANGENDIJKCOMP Comparison process according to Langendijk et al. (2002)
%   Usage:        [ out ] = langecomp( in1,in2,s,do,cp )
%
%   Input parameters:
%     in1    :  modified DFT for one target position and both ears
%     in2    :  stored DFT-templates for all positions and both ears
%
%   Output parameters:
%     out    :  probability density function (pdf)
%
%   XXX Description is missing.
%
%   The function accepts the following optional parameters:
%
%     's', s     Standard deviation of transforming Gaussian function.
%                The default value is 2.
%
%     'do',do    Differential order. The default value is 0.
%
%     'std'      Use the 'std' comparison process. This is the default.
%
%     'xcorr'    Use the 'xcorr' comparison process. 
%
%     'bal',bal  balance of left to right channel. The default is 1.
%
%   See also: langendijk
  
  
% AUTHOR: Robert Baumgartner, OEAW Acoustical Research Institute


definput.import={'langendijkcomp'};
[flags,kv]=ltfatarghelper({'do','s','bal'},definput,varargin);
  
if flags.do_std
  
  ptemp=zeros(size(in2,2),size(in2,3)); % initialisation
  for ch=1:size(in2,3)
    for ind=1:size(in2,2)
      if kv.do==0
        z=in1(:,ch)-in2(:,ind,ch);
      else
        z=diff(in1(:,ch),kv.do)-diff(in2(:,ind,ch),kv.do);
      end
      sigma=std(z,1);
      ptemp(ind,ch)=normpdf(sigma,0,kv.s);
    end
  end
  if size(in2,3)==1
    p=ptemp;
  else
    p=(kv.bal*ptemp(:,1)+1/kv.bal*ptemp(:,2))/2; % balance
  end
  p=1/sum(p).*p;
  % p=p/max(p); % normalisation
    
end;
  
if flags.do_xcorr

  ptemp=zeros(size(in2,2),size(in2,3)); % initialisation
  for ch=1:size(in2,3)
    for ind=1:size(in2,2)
      if kv.do==0
        z=corrcoef(in1(:,ch),in2(:,ind,ch));
      else
        z=corrcoef(diff(in1(:,ch),kv.do),diff(in2(:,ind,ch),kv.do));
      end
      z=z(2);
      ptemp(ind,ch)=(z+1)/2;
    end
  end
  if size(in2,3)==1
    p=ptemp;
  else
    p=(kv.bal*ptemp(:,1)+1/kv.bal*ptemp(:,2))/2; % balance
  end
  
end

out=p;


