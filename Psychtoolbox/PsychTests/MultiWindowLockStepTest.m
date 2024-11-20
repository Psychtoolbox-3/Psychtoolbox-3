function MultiWindowLockStepTest(nrwins, separateScreens)
% MultiWindowLockStepTest([nrwins=2][, separateScreens=0]);
%
% Test if and how many parallel asynchronous window flips
% Psychtoolbox can handle on multiple onscreen windows.
%
% This test exercises the asynchronous flip scheduling and
% timestamping, aka the Screen('AsyncFlipBegin'); and
% Screen('AsyncFlipCheckEnd'); functions for parallel scheduling
% of multiple independent bufferswaps on multiple open onscreen
% windows.
%
% 'nrwins' selects how many onscreen windows to drive in parallel.
% If set to [], 2 windows will be used. The i'th window is flipped
% at a target rate of one flip every i'th video refresh interval.
%
% Timestamps are collected for all flips and in the end one plot
% per window shows how well the flips met their requested deadlines.
%
% If 'nrwins' is omitted, a simple test with only two windows is
% conducted. One window flips after 10 seconds, the other after 20
% seconds.
%
% If the optional parameter 'separateScreens' is set to a non-zero
% setting, then the different onscreen windows are not opened
% on the same display screen 0, but each one on a separate screen.
% Performance can differ significantly between one screen and
% multi-screen, as the operating system will use very different
% underlying algorithms and methods to handle flips on windows,
% depending on their placement on different displays or not. Each
% operating system will also behave differently.
%

% History:
% 07.09.2011  mk  Written.

AssertOpenGL;

if nargin == 0
  w1 = Screen('Openwindow', 0, [255 0 0], [0 0 200 200]);
  w2 = Screen('Openwindow', 0, [0 255 0], [210 0 410 200]);

  DrawFormattedText(w1, 'I swap in 20 seconds.', 'center', 'center');
  DrawFormattedText(w2, 'I swap in 10 seconds.', 'center', 'center');

  tic;
  t1 = Screen('Flip', w1);
  t2 = Screen('Flip', w2);

  Screen('AsyncFlipBegin', w1, t1 + 20, [], 0);
  w1Engaged = 1
  toc
  % Need to force-enable OpenGL context of window w2, or error checks in
  % AsyncFlipBegin will fail atm., as we don't have an OpenGL context bound
  % after initial AsyncFlipBegin. Normally any kind of other Screen()
  % drawing command would prevent this problem...
  Screen('GetWindowInfo', w2);
  Screen('AsyncFlipBegin', w2, t2 + 10, [], 0);
  w2Engaged = 1
  toc

  td2 = Screen('AsyncFlipEnd', w2) - t2
  td1 = Screen('AsyncFlipEnd', w1) - t1
  allfinished = 1
  toc

  KbStrokeWait;

  sca;
  return;
end

% Assign defaults:
if isempty(nrwins) || nrwins < 0 || nrwins > 20
  nrwins = 2;
end

if nargin < 2 || isempty(separateScreens)
  separateScreens = 0;
end

% Init timestamps and counts:
c = zeros(nrwins, 1);
ifi = zeros(nrwins, 1);
t = nan(nrwins, 5000);

screens = Screen('Screens');
if separateScreens
  % One screen per fullscreen window:
  if length(screens) < nrwins
    error('You asked to display each window on a separate screen, but there are less screens than windows!');
  end

  for i = 1:nrwins
    w(i) = Screen('Openwindow', screens(i + double(IsWin)), 0, []);
    ifi(i) = Screen('GetFlipInterval', w(i));
  end
else
  % All tiny windows on first screen:
  for i = 1:nrwins
    w(i) = Screen('Openwindow', 0, 0, [(i-1)*50, 0, i*50, 50]);
    ifi(i) = Screen('GetFlipInterval', w(i));
  end
end

% Wait for go from user:
KbReleaseWait;

% Get baseline timestamp of vsync:
tvbl = Screen('Flip', w(1));
for i = 1:nrwins
  t(i, 1) = tvbl; 
end

% Run until keypress or 5000 frames elapsed:
while ~KbCheck && (c(1) < 5000)
  % Check all windows:
  for i = nrwins:-1:1
    % Check/Process i'th window:
    if Screen('GetWindowInfo', w(i), 4) == 0
      % w(i) doesn't have an async flip scheduled.
      % Render something into it, schedule one:
      c(i) = c(i) + 1;
      x = mod(c(i), 50);
      Screen('DrawLine', w(i), 255, 0, x, 50, x, 1);
      Screen('AsyncFlipBegin', w(i), t(i, c(i)) + (i - 0.5) * ifi(i));
    else
      % w(i) has an async flip pending. Check if completed:
      tvbl = Screen('AsyncFlipCheckEnd', w(i));
      if tvbl > 0
        % Completed - and finalized. Store its timestamp:
        t(i, c(i) + 1) = tvbl;
      end
    end
  end
  WaitSecs('YieldSecs', 0.002);
end

% Close displays:
sca;

% Do the stats and plotting of fancy graphs:
close all;
for i=1:nrwins
  figure;
  tvec = t(i, 1:c(i)-1);
  tvec = diff(tvec) / ifi(i);
  plot(tvec);
  title(sprintf('Should be %i frames delta:', i));
  fprintf('Window %i should have mean %i, has %f with range %f.\n', i, i, mean(tvec), psychrange(tvec));
end

return;
