function ErrorCatchDemo
%
% Show how to use try and catch to clean up after an error.
%
% 12/19/03  dhb  Wrote it based on suggestion from Scott French.
% 10/10/06  rhh  Edited it.
% 10/11/06	dhb   Turn off warnings, use max screen.

try
    % Enclose all your real code between try and catch.
    
    % Removes the blue screen flash and minimize extraneous warnings.
	oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);

    % Find out how many screens and use largest screen number.
    whichScreen = max(Screen('Screens'));
    
    Screen('OpenWindow', whichScreen);
    error('Oops!');
    
    % Clean up, although in this case we never get here.
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

catch
    % If an error occurs, the catch statements executed.  We restore as
    % best we can and then rethrow the error so user can see what it was.
    
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    
    fprintf('We''ve hit an error.\n');
    psychrethrow(psychlasterror);
    fprintf('This last text never prints.\n');
    
end
