
tests_todo={ 'adaptloop', 'outermiddle', 'coefgtdrnl','dbsplsafety'};


total_tests_failed=0;
list_of_failed_tests={};

for ii=1:length(tests_todo)
  test_failed=feval(['test_',tests_todo{ii}]);
  total_tests_failed=total_tests_failed+test_failed;
  if test_failed>0
    list_of_failed_tests{end+1}=['test_',tests_todo{ii}];
  end;
end;

amt_disp(' ');
if total_tests_failed==0
  amt_disp('ALL TESTS PASSED');
else
  s=sprintf('%i TESTS FAILED',total_tests_failed);
  amt_disp(s);
  amt_disp('The following test scripts contained failed tests');
  for ii=1:length(list_of_failed_tests)
    amt_disp(['   ',list_of_failed_tests{ii}]);
  end;
end;


