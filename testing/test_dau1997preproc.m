insig=greasy;
fs=16000;

[outsig_ref, fc_ref, mfc_ref] = ref_dau1997_preproc(insig, fs);

% The reference must use a 'basef', so also do this here.
[outsig, fc, mfc] = dau1997_preproc(insig, fs, 'basef', 1000);

norm(fc-fc_ref)
norm(mfc-mfc_ref)

for nfc=1:length(fc)
  for nmfc=1:size(outsig{nfc},2)
    res=norm(outsig{nfc}(:,nmfc)-outsig_ref(:,nfc,nmfc))/norm(outsig{nfc}(:,nmfc));
    
    amt_disp(sprintf('Freq: %f MFreq: %f res: %f',fc(nfc),mfc(nmfc),res));
  end;  
end;


