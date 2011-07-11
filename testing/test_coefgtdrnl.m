function test_failed=test_coefgtdrnl
%TEST_COEFGTDRNL
%
%  Test how to replace the coefGtDRNL function by gammatone.

disp(' ===============  TEST_COEFGTDRNL ================');

test_failed=0;

fc=2000;
fs=10000;
n=4;
bw=200;

[b1,a1]=coefGtDRNL(fc,bw,n,fs);
[b2,a2]=gammatone(fc,fs,n,bw/audfiltbw(fc),'classic');

res=norm(b1-b2);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GTDRNL B %f %s\n', res,fail);

res=norm(a1-a2);
[test_failed,fail]=ltfatdiditfail(res,test_failed);
fprintf('GTDRNL A %f %s\n', res,fail);
