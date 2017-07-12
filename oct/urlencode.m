function utf=urlencode(str)
% This function emulates the MATLAB function urlencode by encoding 
% the string STR to UTF-8. All special characters are encoded to 
% %HH (ASCII in hex), except: `.` `_` `*` and `-`. Further, ` ` is encoded to `+`.

% 12.07.2017, Piotr Majdak
  utf = '';
  for ii = 1:length(str),
    if isalnum(str(ii)) || str(ii)=='.' || str(ii)=='_' || str(ii)=='*' || str(ii)=='-' 
      utf(end+1) = str(ii);
    elseif str(ii)==' '
      utf(end+1) = '+';
    else
      utf=[utf,'%',dec2hex(str(ii)+0)];
    end; 	
  end  
