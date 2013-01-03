function PerceptualVBLSyncTest(screen, stereomode, fullscreen, doublebuffer, maxduration, vblSync, testdualheadsync)
% PerceptualVBLSyncTest(screen, stereomode, fullscreen, doublebuffer, maxduration, vblSync, testdualheadsync)
%
% Perceptual synchronization test for synchronization of Screen('Flip') and
% Screen('WaitBlanking') to the vertical retrace.
%
% Arguments:
% 'screen' Either a single screen handle, or none (in which case the
% display with the maximum id will be used), or a vector of two handles in
% stereomode 10, e.g., [0 1] if you want to output to screens 0 and 1. You
% can also pass a vector of two screens when stereomode is not set to 10.
% In this case two separate (non-stereo) onscreen windows will be opened on
% both displays and they will get flipped in multiflip mode 2. That means
% that the first display (first element of 'screen') is synced to VBL, but
% the 2nd one is synced to bufferswaps of the first one. This is a
% straightforward test to check if two displays of a stereosetup run with a
% synchronized retrace cycle (good!) or if they are phase-shifted or
% drifting against each other (not good!).
%
% 'stereomode' Which stereomode to use? Defaults to zero, ie. no stereo.
%
% 'fullscreen' Fullscreen presentation? Defaults to 1 ie. yes. In
% non-fullscreen mode, no proper synchronization of bufferswaps can be
% expected.
%
% 'doublebuffer' Single- or double-buffering (1). Defaults to 1. In single
% buffer mode there is no sync to retrace, so this is a good way to
% simulate the tearing artifacts that would happen on sync failure, just to
% get an impression.
%
% 'maxduration' Maximum runtime of test: Runs until keypress or maxduration
% seconds have elapsed (Default is 10 seconds).
%
% 'vblSync' If 1, synchronize bufferswaps to vertical retrace of monitor,
% otherwise (setting 0) swap immediately without sync, ie., usually with tearing.
%
% 'testdualheadsync' If 1, and 'vblSync' is zero, manually wait until the video
% scanout position reaches half the height of the display, then swap. If this
% is done on a multi-display setup and the video scanout cycles of all the
% participating displays are properly synchronized, you should see a "static"
% crack line roughly at half the height of the display, maybe a bit lower. If
% you see a wandering crack line, at least on some displays, or you see vertical
% offsets of the position of the crack line between displays, then the displays
% are not properly synchronized, ie., not suitable for artifact free binocular
% stimulation. Caveat: This logic has been developed and tested specifically
% for testing on Linux with a single X-Screen spanning multiple displays. It may
% or may not be suitable to assess other operating systems or display configurations.
%
% After starting this test, you should see a flickering greyish background
% that flickers in a homogenous way - without cracks or weird moving patterns
% in the flickering area. If you see an imhogenous flicker, this means that
% synchronization of stimulus onset to the vertical retrace doesn't work due
% to some serious bug or limitation of your graphics hardware or its driver.
% If you don't know what this means, you can test this script with parameter
% doublebuffer == 0 to artificially create a synchronization failure.
%
% On many systems you should also see some emerging pattern of yellow horizontal lines.
% These lines should be tightly concentrated/clustered in the topmost area of
% the screen. Lots of yellow lines in the middle or bottom area or even
% randomly distributed lines indicate some bug in the driver of your graphics
% hardware. This is a common problem of all ATI graphics adapters on MacOS-X
% versions earlier than OS-X 10.4.3 when running a dual-display setup...
%
% A second reason for distributed yellow lines could be bad timing on your
% machine, e.g., due to background activity by virus scanners or the Spotlight
% indexing service on OS-X. Turn these off for conducting your studies!
%

% History:
% 01/28/06 mk Written. Replaces the built-in flickertest of Screen('OpenWindow')
%             on multi-display setups. That test wasn't well received :(

% Check for presence of OpenGL PTB.
AssertOpenGL;

if nargin < 4
    doublebuffer = [];
end

if isempty(doublebuffer)
   % Use double-buffered windows by default: Single-buffered ones can't sync to
   % retrace and are discouraged anyway. Setting doublebuffer=0 is an easy way
   % to reproduce the visual pattern created by a complete sync-failure though.
   doublebuffer=1;
end;
doublebuffer=doublebuffer+1;

if nargin < 2
    stereomode = [];
end

if isempty(stereomode)
   % Use non-stereo display by default. 
   stereomode=0;
end;

if nargin < 3
    fullscreen = [];
end

if isempty(fullscreen)
   fullscreen=1;
end;

if nargin < 1
    screen = [];
end

if isempty(screen)
    if stereomode == 10
        screen(1) = max(Screen('Screens')) - 1;
        screen(2) = max(Screen('Screens'));
        if screen(1)<0
            error('Stereomode 10 only works on setups with two attached displays!');
        end
    else
        screen=max(Screen('Screens'));
    end
end;

if nargin < 5
    maxduration = [];
end

if isempty(maxduration)
    maxduration = 10;
end

if nargin < 6
    vblSync = [];
end

if isempty(vblSync)
    vblSync = 1;
end

if nargin < 7
    testdualheadsync = [];
end

if isempty(testdualheadsync)
    testdualheadsync = 0;
end

thickness = (1-vblSync) * 4 + 1;

try
    if fullscreen
        rect1=[];
        rect2=[];
    else
        rect1=InsetRect(Screen('GlobalRect', screen(1)), 1, 0);
        if length(screen)>1
            rect2=InsetRect(Screen('GlobalRect', screen(2)), 1, 0);
        end
    end;
   
   help PerceptualVBLSyncTest;
   fprintf('Press ENTER key to start the test. The test will stop after 10 seconds\n');
   fprintf('or any keypress...\n');

   KbStrokeWait;
   
   if stereomode~=10
       % Standard case:
       [win , winRect]=Screen('OpenWindow', screen(1), 0, rect1, [], doublebuffer, stereomode);
       if length(screen)>1
           win2 = Screen('OpenWindow', screen(2), 0, rect2, [], doublebuffer, stereomode);
       end
   else
       % Special case for dual-window stereo:

       % Setup master window:
       [win , winRect]=Screen('OpenWindow', screen(1), 0, rect1, [], doublebuffer, stereomode);
       % Setup slave window:
       Screen('OpenWindow', screen(2), 0, rect2, [], doublebuffer, stereomode);       
   end
   
   flickerRect = InsetRect(winRect, 100, 0);
   color = 0;
   deadline = GetSecs + maxduration;
   beampos=0;
   
   ifi = Screen('GetFlipInterval', win);
   
   VBLTimestamp = Screen('Flip', win, 0, 2);
   
   while (~KbCheck) && (GetSecs < deadline)
      % Draw left eye view (if stereo enabled):
      Screen('SelectStereoDrawBuffer', win, 0);
      % Draw alternating black/white rectangle:
      Screen('FillRect', win, color, flickerRect);
      % If beamposition is available (on OS-X), visualize it via yellow horizontal line:
      if (beampos>=0), Screen('DrawLine', win, [255 255 0], 0, beampos, winRect(3), beampos, thickness); end;
      % Same for right-eye view...
      Screen('SelectStereoDrawBuffer', win, 1);
      Screen('FillRect', win, color, flickerRect);
      if (beampos>=0), Screen('DrawLine', win, [255 255 0], 0, beampos, winRect(3), beampos, thickness); end;
      
      if stereomode == 0 && length(screen)>1
          Screen('FillRect', win2, color, flickerRect);
          Screen('DrawingFinished', win2, 0, 2);
          Screen('DrawingFinished', win, 0, 2);
          multiflip = 2;
      else
          multiflip = 0;
      end
      
      % Alternate drawing color from white -> black, or black -> white
      color=255 - color;
      
      if doublebuffer>1
          if vblSync
              % Flip buffer on next vertical retrace, query rasterbeam position on flip, if available:
              [VBLTimestamp, StimulusOnsetTime, FlipTimestamp, Missed, beampos] = Screen('Flip', win, VBLTimestamp + ifi/2, 2, [], multiflip);
          else
              if testdualheadsync == 1
                  Screen('DrawingFinished', win, 0, 1);
                  beampos = -1000;
                  while abs(beampos - winRect(4)/2) > 5
                      beampos = Screen('GetWindowInfo', win, 1);
                  end
              end

              % Flip immediately without sync to vertical retrace, do clear
              % backbuffer after flip for visualization purpose:
              VBLTimestamp = Screen('Flip', win, VBLTimestamp + ifi/2, 0, 2, multiflip);
              % Above flip won't return a 'beampos' in non-VSYNC'ed mode,
              % so we query it manually:
              beampos = Screen('GetWindowInfo', win, 1);
              
              % Throttle a little bit for visualization purpose:
              WaitSecs('YieldSecs', 0.005);
          end
      else
          % Just wait a bit in non-buffered case:
          pause(0.001);
      end;
   end;
   
   Screen('CloseAll');
   return;   
catch
   Screen('CloseAll');   
end;
