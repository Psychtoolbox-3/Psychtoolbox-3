function DrawTextFontSwitchSpeedTest
% DrawTextFontSwitchSpeedTest - Test font settings switching speed.
%
% This is a benchmark of the speed with which Screen('DrawText') et al.
% can handle switching between different font settings between text
% drawing. Specifically it tests switching between different text sizes, text
% styles and font families. It also performs non-switched drawing to establish
% a baseline.
%
% After running for a thousand repetitions for each conditon, it prints a
% summary of drawing speed and then exits.
%
% Background:
%
% In an ideal world, changing text settings should not slow down text drawing.
% In the real world, each time you change a text setting, Screen() needs to
% recompute various things, e.g., it needs to load new font data from the
% filesystem, rasterize glyphs into little images, convert them into OpenGL
% textures and upload them into the graphics card. This takes time and slows
% down text drawing after each change of settings.
%
% The Linux text renderer can cache the precomputed data for multiple different
% sets of text settings. This means that if you alternate repeatedly between a limited
% number of different text settings, e.g., drawing text strings of different font, size or
% style within one stimulus frame, then this will be only slow on the first stimulus frame.
% Linux ptb will remember the computed font data in a cache and recycle it in successive
% repetitions of those text settings. Therefore any repetition of a previously used text
% setting is "free" and fast. The number of cacheable settings is limited, as it consumes
% memory, currently to about 40 different settings. If you exceed that number, Screen
% will expunge the least used settings and be intelligent about what to recycle and what
% not. If you use non-repetitive settings then of course no caching in the world will help
% you. If no caching can be used - or at first use of a setting - subsequent text drawing will
% be about 10 - 15 times slower.
%
% On other operating systems, nothing will be cached or recycled, resulting in generally
% much slower text drawing speeds, even with constant text settings.
%

% History:
% 22-Sep-2013  mk    Written.

PsychDefaultSetup(1);

screenid = max(Screen('Screens'));
w = Screen('OpenWindow', screenid);
n = 1000;

% Test size changes:
textSize = [48, 49];
Screen('TextSize', w, textSize(1));
Screen('DrawText', w, 'Hello World!', 0, 0);

for switchit = 0:1
  t1 = GetSecs;
  for i=1:n
    if switchit
      Screen('TextSize', w, textSize(mod(i,2) + 1));
    else
      Screen('TextSize', w, textSize(1));
    end
    
    Screen('DrawText', w, 'Hello World!', 0, 0);
  end
  t2 = GetSecs;
  Screen('Flip', w);

  if switchit
    fprintf('DrawText speed for size-change: %f msecs/draw.\n', (t2 - t1) * 1000 / n);
  else
    fprintf('DrawText speed for constant size: %f msecs/draw.\n', (t2 - t1) * 1000 / n);
  end
end

% Test style changes:
Screen('TextStyle', w, 0);

for switchit = 0:1
  t1 = GetSecs;
  for i=1:n
    if switchit
      Screen('TextStyle', w, mod(i,2));
    else
      Screen('TextStyle', w, 0);
    end
    
    Screen('DrawText', w, 'Hello World!', 0, 0);
  end
  t2 = GetSecs;
  Screen('Flip', w);

  if switchit
    fprintf('DrawText speed for style-change: %f msecs/draw.\n', (t2 - t1) * 1000 / n);
  else
    fprintf('DrawText speed for constant style: %f msecs/draw.\n', (t2 - t1) * 1000 / n);
  end
end

% Test font changes:
textFont = {'Courier New', 'Times'};
Screen('TextFont', w, textFont{1});
Screen('DrawText', w, 'Hello World!', 0, 0);

for switchit = 0:1
  t1 = GetSecs;
  for i=1:n
    if switchit
      Screen('TextFont', w, textFont{mod(i,2) + 1});
    else
      Screen('TextFont', w, textFont{1});
    end
    
    Screen('DrawText', w, 'Hello World!', 0, mod(i,2) * 100);
  end
  t2 = GetSecs;
  Screen('Flip', w);

  if switchit
    fprintf('DrawText speed for font-change: %f msecs/draw.\n', (t2 - t1) * 1000 / n);
  else
    fprintf('DrawText speed for constant font: %f msecs/draw.\n', (t2 - t1) * 1000 / n);
  end
end

sca;
end
