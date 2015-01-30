
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

amtdisp(' ');
if total_tests_failed==0
  amtdisp('ALL TESTS PASSED');
else
  s=sprintf('%i TESTS FAILED',total_tests_failed);
  amtdisp(s);
  amtdisp('The following test scripts contained failed tests');
  for ii=1:length(list_of_failed_tests)
    amtdisp(['   ',list_of_failed_tests{ii}]);
  end;
end;


