function FontDemo
% FontDemo
%
% Shows how to get the font you want, when you're not
% sure what's available, and only certain fonts will do.  This
% is a portable solution to the problem of running on various
% computers and Operating Systems that offer various fonts.
%
% In the FontDemo program we call CHAR to convert a one-string list,
% e.g. {'Arial'}, to a string, e.g. 'Arial'.

% 5/23/02	 dgp   Wrote it.
% 10/11/06	 dhb   Turn off warnings, use max screen.
% 10/14/2006 dhb  Save and restore altered prefs, more extensive comments for them.
% 09/22/1013 mk   Make robust against lazy-evaluating text renderes, e.g., Linux.

try

    % Screen is able to do a lot of configuration and performance checks on
    % open, and will print out a fair amount of detailed information when
    % it does.  These commands supress that checking behavior and just let
    % the demo go straight into action.  See ScreenTest for an example of
    % how to do detailed checking.
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    
    % Find out how many screens and use largest screen number.
    whichScreen = max(Screen('Screens'));
   
    % Open a new window.
    window = Screen('OpenWindow', whichScreen); 
    fprintf('New window opened with font initially requested to be ''%s''.\n', Screen('TextFont', window));
    % Perform dummy draw call, to take deferred font setting on some OS'es into account:
    Screen('DrawText', window, ' ', 0, 0);
    fprintf('New window opened with font actual initial setting ''%s''.\n', Screen('TextFont', window));
    
    for font = {'Arial', 'Geneva', 'Helvetica', 'sans-serif', 'DejaVu Sans'}
        % Set font
        Screen('TextFont', window, char(font));     % Need "char" input.

        % Perform dummy draw call, to take deferred font setting on some OS'es into account:
        Screen('DrawText', window, ' ', 0, 0);
        
        % Check if the screen's font has changed to the one we requested.
        match = streq(font, Screen('TextFont', window));
        if match
            fprintf('Successfully set font to ''%s''.\n', Screen('TextFont', window));
            break;
        end
    end
    
    if ~match
        fprintf('The real window font is ''%s''.\n', Screen('TextFont', window));
        error('We couldn''t find any of the fonts we wanted.');
    end
    
    Screen('Close', window);
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
catch
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    psychrethrow(psychlasterror);
end
