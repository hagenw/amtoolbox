function op = do_xcorr(left, right, fs, fc)
iacc = xcorr(left,right,round(fs/(fc*2)),'coeff');
[coherence delay_samp] = max(iacc);
phase = fc*2*pi*delay_samp/fs;
op = [phase coherence];