function SimpleMovieDemo(moviename, windowrect)
% Most simplistic demo on how to play a movie.
%
% SimpleMovieDemo(moviename [, windowrect=[]]);
%
% This bare-bones demo plays a single movie whose name has to be provided -
% including the full filesystem path to the movie - exactly once, then
% exits. This is the most minimalistic way of doing it. For a more complex
% demo see PlayMoviesDemo. The remaining demos show more advanced concepts
% like proper timing etc.
%
% The demo will play our standard DualDiscs.mov movie if the 'moviename' is
% omitted.
%

% History:
% 02/05/2009  Created. (MK)
% 06/17/2013  Cleaned up. (MK)

% Check if Psychtoolbox is properly installed:
AssertOpenGL;

if nargin < 1 || isempty(moviename)
    % No moviename given: Use our default movie:
    moviename = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
end

if nargin < 2 || isempty(windowrect)
    windowrect = [];
end

% Wait until user releases keys on keyboard:
KbReleaseWait;

% Select screen for display of movie:
screenid = max(Screen('Screens'));

try
    % Open 'windowrect' sized window on screen, with black [0] background color:
    win = Screen('OpenWindow', screenid, 0, windowrect);
    
    % Open movie file:
    movie = Screen('OpenMovie', win, moviename);
    
    % Start playback engine:
    Screen('PlayMovie', movie, 1);
    
    % Playback loop: Runs until end of movie or keypress:
    while ~KbCheck
        % Wait for next movie frame, retrieve texture handle to it
        tex = Screen('GetMovieImage', win, movie);
        
        % Valid texture returned? A negative value means end of movie reached:
        if tex<=0
            % We're done, break out of loop:
            break;
        end
        
        % Draw the new texture immediately to screen:
        Screen('DrawTexture', win, tex);
        
        % Update display:
        Screen('Flip', win);
        
        % Release texture:
        Screen('Close', tex);
    end
    
    % Stop playback:
    Screen('PlayMovie', movie, 0);
    
    % Close movie:
    Screen('CloseMovie', movie);
    
    % Close Screen, we're done:
    Screen('CloseAll');
    
catch %#ok<CTCH>
    sca;
    psychrethrow(psychlasterror);
end

return
