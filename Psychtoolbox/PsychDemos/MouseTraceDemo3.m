% MouseTraceDemo3
%
% Allow multiple people to draw simultaneously by use of multiple
% mice or other pointing devices, utilizing Multi-Pointer-X on Linux.
%
% See also: PsychDemos, MouseTraceDemo, GetMouse.
%
% HISTORY
%                       
% 7/30/11  mk       Derived from MouseTraceDemo2.
%

AssertOpenGL;

if ~IsLinux
  error('Sorry, this demo currently only works on Linux.');
end

try
    % Open up a window on the screen and clear it.
    whichScreen = max(Screen('Screens'));
    [theWindow,theRect] = Screen('OpenWindow', whichScreen, 0);

    % Use an offscreen window as drawing canvas:
    woff = Screen('OpenOffscreenWindow', theWindow, 0);

    % Get handles for all virtual pointing devices, aka cursors:
    mice = GetMouseIndices('masterPointer');

    % Move the virtual cursors to the center of the screen
    for mouse = mice
      theX(mouse+1) = theRect(RectRight)/2;
      theY(mouse+1) = theRect(RectBottom)/2;
      SetMouse(theX(mouse+1), theY(mouse+1), whichScreen, mouse);
      col(mouse+1, 1:3) = rand(1,3) * 255;

      % Hide the system-generated cursors. We do this, because only the
      % first mouse cursor is hardware-accelerated, ie., a GPU created
      % hardware cursor. All other cursors are software-cursors, created
      % by the Windowing system. These tend to flicker badly in our use
      % case. Therefore we disable all system cursor images and draw our
      % cursors ourselves for a more beautiful look:
      HideCursor([], mouse);
    end

    % Some instructions, drawn into the drawing canvas:
    Screen(woff,'FillRect',0);
    Screen(woff,'TextSize',24);
    Screen(woff,'DrawText','Drag mice (i.e. hold button down) to draw. Any key for exit.',50,50,255);

    % Wait for release of all keys on all keyboards:
    KbReleaseWait(-3);

    % Stay in redraw loop as long as no key on any keyboard pressed:
    while ~KbCheck(-3)
      % Blit offscreen window with users scribbling into onscreen window:
      Screen('DrawTexture', theWindow, woff);

      % Check all masterpointer mouse/pointing devices:
      for mouse = mice
        [x(mouse+1), y(mouse+1), buttons] = GetMouse(theWindow, mouse);
	if any(buttons(1:3))
	  % Update offscreen window with latest scribbling from user for this 'mouse':
          if (x(mouse+1) ~= theX(mouse+1) || y(mouse+1) ~= theY(mouse+1))
            Screen('DrawDots', woff, [x(mouse+1), y(mouse+1)], 10, col(mouse+1, :));
            theX(mouse+1) = x(mouse+1); theY(mouse+1) = y(mouse+1);
          end
	end

	% Draw a dot to visualize the mouse cursor for this 'mouse:
	Screen('DrawDots', theWindow, [x(mouse+1), y(mouse+1)], 5, col(mouse+1, :));
      end

      % Flip the updated onscreen window:
      Screen('Flip', theWindow);
    end

    % Show master cursors again:
    for mouse = mice
      ShowCursor('Arrow', [], mouse);
    end

    Screen('CloseAll');
catch
    sca;
    psychrethrow(psychlasterror);
end %try..catch..
