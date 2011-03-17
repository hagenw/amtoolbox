function [test_failed,fail]=amtdiditfail(res,test_failed);
%AMTDIDITFAIL  Did a test fail
%
%  [test_fail,fail]=AMTDIDITFAIL(res,test_fail) updates test_fail if
%  res is above threshhold and outputs the word FAIL in the variable
%  fail. Use only in testing scripts.
 
fail='';
if res>10e-10
  fail='FAILED';
  test_failed=test_failed+1;
end;
