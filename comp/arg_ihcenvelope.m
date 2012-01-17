function definput=arg_ihcenvelope(definput)
 
  definput.flags.ihctype={'nodefault','ihc_bernstein','ihc_breebaart','ihc_dau','hilbert', ...
                    'ihc_lindemann','ihc_meddis'};

  definput.keyvals.minlvl=[];
