function op = calc_benefit(target_in,int_in,ears,fs)

nerbs = 1:1:round(erb_rate(fs/2));
effective_SNR = zeros(1,length(nerbs)); fc = zeros(size(nerbs)); bmld_prediction = zeros(size(nerbs));

for n = 1:length(nerbs)
    fc(n) = round(f_of_erb_rate(nerbs(n)));        % get filter cf 
% filter target and interferer separately
    targ_left = gammatone_c(target_in(:,1),fs,fc(n));    
    targ_right = gammatone_c(target_in(:,2),fs,fc(n));   
    int_left = gammatone_c(int_in(:,1),fs,fc(n));       
    int_right = gammatone_c(int_in(:,2),fs,fc(n));
% BMLD
    if ears == 'both'
      int_stats = do_xcorr(int_left,int_right,fs,fc(n)); % cross-correlate
      targ_stats = do_xcorr(targ_left,targ_right,fs,fc(n));    
      bmld_prediction(n) = bmld(int_stats(2),targ_stats(1),int_stats(1),fc(n));
    end
% better-ear SNR
    right_SNR = sum(targ_right.^2) / sum(int_right.^2);
    left_SNR = sum(targ_left.^2) / sum(int_left.^2);
    if ears == 'both'
      SNR = max(left_SNR,right_SNR);
    elseif ears == 'left'
      SNR = left_SNR;
    else
      SNR = right_SNR;
    end
   
% combination
    effective_SNR(n) = 10*log10(SNR);
end

weightings = calc_weightings(fc);
if ears == 'both'
  weighted_bmld = sum(bmld_prediction.*weightings);
else
  weighted_bmld = 0;
end
weighted_SNR = sum(effective_SNR.*weightings);
benefit = weighted_SNR + weighted_bmld;

op=[benefit weighted_SNR weighted_bmld];

