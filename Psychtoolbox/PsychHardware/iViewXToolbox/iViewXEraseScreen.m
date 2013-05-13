function iViewXEraseScreen(ivx)

% OS X version
if ~isempty(ivx.window)
Screen('FillRect',ivx.window, ivx.backgroundColour);
Screen('Flip',ivx.window);
else
    fprintf([mfilename ': cannot erase screen, no window defined\n']);
end
