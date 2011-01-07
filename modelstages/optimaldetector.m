function [Y,mue1,mue2,mue3] = optimaldetector(repdiff,repdet)
%OPTIMALDETECTOR  Generic optimal detector for the CASP and Breebart models
%
  
corr_mue = repdiff .* repdet;
dims = size(corr_mue);
optfactor = sqrt(prod(dims));

mue = mean(corr_mue,1);
mue1 = squeeze(mue);

if dims(2) >= 2
  mue=mean(mue,2);
  mue2 = squeeze(mue);
end

if (ndims(corr_mue)>2) && (dims(3) >= 2)
  mue=mean(mue,3);
  mue3 = mue;
end

mue = mue * optfactor;
Y = mue;
