function test_failed=test_outermiddle
%TEST_DAU96  Test outermiddle
%
%  Test outermiddel by comparison with reference implementation.

disp(' ===============  TEST_OUTERMIDDLE ================');
  
test_failed=0;
  
fs=22050;

L=10*fs;

f=randn(L,1);

ref_HeadphoneFilter(fs);
ref_MiddleEarFilter(fs);

fr=OuterMiddleFilter(f);

outer_ear_fir_coeff=headphonefilter(fs);
mid_ear_fir_coeff=middleearfilter(fs,'jepsenmiddleear');

ff = fftfilt(outer_ear_fir_coeff,f);
ff = fftfilt(mid_ear_fir_coeff,ff);

res=norm(ff(:)-fr(:));

[test_failed,fail]=ltfatdiditfail(res,test_failed);

disp(sprintf('OUTER+MIDDLE %0.5g %s',res,fail));

