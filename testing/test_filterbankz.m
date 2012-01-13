nsigs=3;

blocksize=4096;

nblocks=5;

fs=16000;

hopsize=2;

siglen=nblocks*blocksize+10*hopsize;



[b,a] = gammatone(erbspacebw(0,fs/2),fs,'complex');

fb=filterbank_init(b,a,nsigs,hopsize);

insig  = randn(siglen,nsigs);
outsig = zeros(siglen/hopsize,fb.nchannels,nsigs);

for ii=0:nblocks
  [outsig_block,fb]=filterbank_block(insig(1+ii*blocksize:min((ii+1)* ...
                                           blocksize,siglen),:),fb);
  
  outsig(fb.outstart:fb.outend,:,:)=outsig_block;
    

end;

outsig_ref=ufilterbankz(b,a,insig,hopsize);

norm(outsig(:)-outsig_ref(:))



%OLDFORMAT
