function TextFontTest

% Test setting the text font.

try;
    w=Screen('OpenWindow',0);
%     Screen(w,'TextFont','Chicago');
    Screen('TextFont',w,'Chicago');
    clear screen
catch;
    clear screen
%     rethrow(lasterror);
end
