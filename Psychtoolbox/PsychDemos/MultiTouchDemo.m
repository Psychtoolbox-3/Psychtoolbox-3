function MultiTouchDemo(dev, verbose)
% MultiTouchDemo([dev][, verbose=0]) - A advanced demo for multi-touch touchscreens.
%
% Run it. Pressing the ESCape key will stop it.
%
% Touch the screen and watch the nice colorful happy blobs
% sprinkle to life :)
%
% The demo will try to use the first available touchscreen, or if
% there isn't any touchscreen, the first available touchpad. You
% can also select a specific touch device by passing in its 'dev'
% device handle. Use of touchpads usually needs special configuration.
% See "help TouchInput" for more info.
%
% If you set 'verbose' to 1, then all the available info about each
% touch point will be displayed close to the touch point. As drawing
% so much text is very slow, the demo will only update touchpoints
% very slowly!
%
% The demo not only tracks and display touches and their location, as
% well as their timing. It also uses info about the shape of the contact,
% visualizing this as aspect ratio of the drawn rectangle. It tries to
% get info about the contacts orientation and draws accordingly. It tries
% to get (or derive) info about touch pressure and modulates the brightness
% of the rectangle accordingly. One property not used, but available would
% be distance for fingers or tools hovering over a touch surface, if that
% surface provides that info.
%
% This demo currently only works on Linux + X11 display system,
% not on Linux + Wayland, not on other operating systems.
%
% For background info on capabilities and setup see "help TouchInput".
%

% History:
% 01-Oct-2017 mk  Written.
% 03-Aug-2018 mk  Only exit demo on ESC key, not on all keys. Doc fixes.

  % Setup useful PTB defaults:
  PsychDefaultSetup(2);

  if nargin < 1
    dev = [];
  end

  if nargin < 2 || isempty(verbose)
    verbose = 0;
  end

  % If no user-specified 'dev' was given, try to auto-select:
  if isempty(dev)
    % Get first touchscreen:
    dev = min(GetTouchDeviceIndices([], 1));
  end

  if isempty(dev)
    % Get first touchpad:
    dev = min(GetTouchDeviceIndices([], 0));
  end

  if isempty(dev) || ~ismember(dev, GetTouchDeviceIndices)
    fprintf('No touch input device found, or invalid dev given. Bye.\n');
    return;
  else
    fprintf('Touch device properties:\n');
    info = GetTouchDeviceInfo(dev);
    disp(info);
  end

  % Open a default onscreen window with black background color and
  % 0-1 color range:
  [w, rect] = PsychImaging('OpenWindow', 0, 0)
  baseSize = RectWidth(rect) / 20;

  % No place for you, little mouse cursor:
  HideCursor(w);

  % Create a 5 pixel texture, just to define a shape we can easily scale and rotate:
  finger = Screen('MakeTexture', w, [1; 0.6; 0.4; 0.4; 0.4]);

  try
    % Create and start touch queue for window and device:
    TouchQueueCreate(w, dev);
    TouchQueueStart(dev);

    % Wait for the go!
    KbReleaseWait;

    % blobcol tracks active touch points - and dying ones:
    blobcol = {};
    buttonstate = 0;
    colmap = [ 1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 1, 0, 1; 0, 1, 1; 1, 1, 1];

    % Only ESCape allows to exit the demo:
    RestrictKeysForKbCheck(KbName('ESCAPE'));

    % Main loop: Run until keypress:
    while ~KbCheck
      % Process all currently pending touch events:
      while TouchEventAvail(dev)
        % Process next touch event 'evt':
        evt = TouchEventGet(dev, w);

        % Touch blob id - Unique in the session at least as
        % long as the finger stays on the screen:
        id = evt.Keycode;

        if evt.Type == 0
          % Not a touch point, but a button press or release on a
          % physical (or emulated) button associated with the touch device:
          buttonstate = evt.Pressed;
          continue;
        end

        if evt.Type == 1
          % Not really a touch point, but movement of the
          % simulated mouse cursor, driven by the primary
          % touch-point:
          Screen('DrawDots', w, [evt.MappedX; evt.MappedY], baseSize, [1,1,1], [], 1, 1);
          continue;
        end

        if evt.Type == 2
          % New touch point -> New blob!
          blobcol{id}.col = colmap(mod(id, 7) + 1, :);
          blobcol{id}.mul = 1.0;
          blobcol{id}.x = evt.MappedX;
          blobcol{id}.y = evt.MappedY;
          blobcol{id}.t = evt.Time;
          blobcol{id}.dt = 0;
        end

        if evt.Type == 3
          % Moving touch point -> Moving blob!
          blobcol{id}.x = evt.MappedX;
          blobcol{id}.y = evt.MappedY;
          blobcol{id}.dt = ceil((evt.Time - blobcol{id}.t) * 1000);
          blobcol{id}.t = evt.Time;
        end

        if evt.Type == 4
          % Touch released - finger taken off the screen -> Dying blob!
          blobcol{id}.mul = 0.999;
          blobcol{id}.x = evt.MappedX;
          blobcol{id}.y = evt.MappedY;
        end

        if evt.Type == 5
          % Lost touch data for some reason:
          % Flush screen red for one video refresh cycle.
          fprintf('Ooops - Sequence data loss!\n');
          Screen('FillRect', w, [1 0 0]);
          Screen('Flip', w);
          Screen('FillRect', w, 0);
          continue;
        end

        if ismember(evt.Type, [2,3,4])
          evt = GetTouchValuators(evt, info);
          if isfield(evt, 'TouchMajor') && isfield(evt, 'TouchMinor')
            % Shape of primary touch ellipse known:
            if evt.TouchMajor > 0 && evt.TouchMinor > 0
              aspect = evt.TouchMajor / evt.TouchMinor;
            else
              aspect = 1;
            end

            blobcol{id}.rect = [0, 0, baseSize, baseSize * aspect];
            if isfield(evt, 'WidthMajor') && isfield(evt, 'WidthMinor')
              % Shape of approach ellipse known: Can use this to approximate
              % pressure if pressure doesn't get reported:
              blobcol{id}.pressure = evt.TouchMajor / evt.WidthMajor;
            else
              blobcol{id}.pressure = 1;
            end
          else
            % Shape unknown, so assume unit shape:
            blobcol{id}.rect = [0 0  baseSize baseSize];
            blobcol{id}.pressure = 1;
          end

          if isfield(evt, 'Orientation')
            blobcol{id}.orientation = evt.Orientation;
          else
            blobcol{id}.orientation = 0;
          end

          if isfield(evt, 'Pressure')
            blobcol{id}.pressure = evt.Pressure;
          end
        end

        if verbose
          blobcol{id}.text = disp(evt);
        end
      end

      % Now that all touches for this iteration are processed, repaint screen
      % with all live blobs at their new positions, and fade out the dying/orphaned
      % blobs:
      for i=1:length(blobcol)
        if ~isempty(blobcol{i}) && blobcol{i}.mul > 0.1
          % Draw it: .mul defines size of the blob:
          Screen('DrawTexture', w, finger, [], CenterRectOnPointd(blobcol{i}.rect * blobcol{i}.mul, blobcol{i}.x, blobcol{i}.y), ...
                 blobcol{i}.orientation, [], [], blobcol{i}.col * blobcol{i}.pressure);
          if blobcol{i}.mul < 1
            % An orphaned blob with no finger touching anymore, so slowly fade it out:
            blobcol{i}.mul = blobcol{i}.mul * 0.95;
          else
            % An active touch. Print its unique touch id and dT timestamp delta between updates in msecs:
            Screen('DrawText', w, num2str(i), blobcol{i}.x, blobcol{i}.y, [1 1 0]);
            Screen('DrawText', w, num2str(blobcol{i}.dt), blobcol{i}.x, blobcol{i}.y - 100, [1 1 0]);
            if verbose
              DrawFormattedText(w, blobcol{i}.text, blobcol{i}.x, blobcol{i}.y + 30, [1 1 0]);
            end
          end
        else
          % Below threshold: Kill the blob:
          blobcol{i} = [];
        end
      end

      if buttonstate
        Screen('FrameRect', w, [1, 1, 0], [], 5);
      end

      % Done repainting - Show it:
      Screen('Flip', w);

      % Next touch processing -> redraw -> flip cycle:
    end

    TouchQueueStop(dev);
    TouchQueueRelease(dev);
    RestrictKeysForKbCheck([]);
    ShowCursor(w);
    sca;
  catch
    TouchQueueRelease(dev);
    RestrictKeysForKbCheck([]);
    sca;
    psychrethrow(psychlasterror);
  end
end
