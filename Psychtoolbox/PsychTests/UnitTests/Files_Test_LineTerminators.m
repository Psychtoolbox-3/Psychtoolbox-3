function success = Files_Test_LineTerminators

success = true;

% create some input
a = sprintf('1\n2');        % unix/linux/OSX
b = sprintf('1\r\n2');      % windows
c = sprintf('1\r2');        % macintosh (OS9 or earlier)
d = sprintf('1\r\n\r2');    % windows, followed by mac os9

try
    % test: ReplaceLineTerminators
    success = tester(ReplaceLineTerminators(a,'unix'),sprintf('1\n2'),'ReplaceLineTerminators(a,''unix'')') && success;
    success = tester(ReplaceLineTerminators(a,'windows'),sprintf('1\r\n2'),'ReplaceLineTerminators(a,''windows'')') && success;
    success = tester(ReplaceLineTerminators(b,'unix'),sprintf('1\n2'),'ReplaceLineTerminators(b,''unix'')') && success;
    success = tester(ReplaceLineTerminators(b,'os9'),sprintf('1\r2'),'ReplaceLineTerminators(b,''os9'')') && success;
    success = tester(ReplaceLineTerminators(c,'unix'),sprintf('1\n2'),'ReplaceLineTerminators(c,''unix'')') && success;
    % this detects one windows terminator and one os9, both get replaced
    % with \n, unix terminator.
    success = tester(ReplaceLineTerminators(d,'unix'),sprintf('1\n\n2'),'ReplaceLineTerminators(d,''unix'')') && success;
    % this detects one windows terminator and one os9, both get replaced
    % with \r\n, windows terminator.
    success = tester(ReplaceLineTerminators(d,'windows'),sprintf('1\r\n\r\n2'),'ReplaceLineTerminators(d,''winsdows'')') && success;
    
    
    % test: BreakLines
    success = tester(BreakLines(a),{'1' '2'},'BreakLines(a)') && success;
    success = tester(BreakLines(b),{'1' '2'},'BreakLines(b)') && success;
    success = tester(BreakLines(c),{'1' '2'},'BreakLines(c)') && success;
    % two terminators in a row, so empty line->cell in between
    success = tester(BreakLines(d),{'1' char(zeros(1,0)) '2'},'BreakLines(d)') && success;
    
catch
    success = false;
    fprintf('error ocurred: "%s"\n',lasterr);
end


function [success] = tester(out,check,funcall)

if ~isequal(out,check)
    success = false;
    fprintf('%s failed\n',funcall);
else
    success = true;
end
