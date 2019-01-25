function i = findnum(str)
% Find values in a string which can be converted to numeric without
% returning an error.
%
% "str" is the string in which to find numbers
%
% "i" is a logical matrix with the indices of numbers in the string as
% true.

i = find(... % Find indices at which str is equal to...
    str == '0' | ... % ...0
    str == '1' | ... % ...1
    str == '2' | ... % ...2
    str == '3' | ... % ...3
    str == '4' | ... % ...4
    str == '5' | ... % ...5
    str == '6' | ... % ...6
    str == '7' | ... % ...7
    str == '8' | ... % ...8
    str == '9' | ... % ...9
    str == '.' | ... % ...a decimal point
    str == '-'   ... % ...a minus sign
    );
end