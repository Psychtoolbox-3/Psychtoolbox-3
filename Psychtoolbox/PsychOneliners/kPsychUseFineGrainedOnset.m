function rc = kPsychUseFineGrainedOnset
% kPsychUseFineGrainedOnset -- Try to create onscreen windows with fine-grained onset scheduling.
%
% This flag can be passed to the optional 'specialFlags' parameter of
% Screen('OpenWindow', ...) or PsychImaging('OpenWindow', ...).
%
% If this flag is set and Psychtoolbox is running on a suitable combination of
% operating system + graphics/display drivers + graphics card + display connector
% + display device and the onscreen window is a non-transparent, topmost,
% unoccluded, decorationless (non-GUI) fullscreen window and everything is set
% up correctly, then a stimulus onset scheduling method will be used that tries
% to achieve that stimulus onsets requested via Screen('Flip', windowPtr, tWhen)
% or Screen('AsyncFlipBegin', windowPtr, tWhen) will no longer be locked in time
% to the start of a fixed refresh rate cycle. In other words, the display will
% no longer be operated at a fixed refresh rate like 60 Hz and true stimulus onset
% will no longer lock onto the start of the next video refresh cycle with a time
% t >= tWhen. Instead the software will try to schedule onset at exactly the
% requested time tWhen, or as closely as possible. This would allow to present
% stimuli at a millisecond granularity target time, instead of a time that is an
% integral multiple of the video refresh rate, e.g., 16.6 msecs increments on a
% 60 Hz display. This would allow more flexible interstimulus intervals or variable
% redraw rates other than, e.g., 60 Hz, 30 Hz, 15 Hz, ...
%
% The technique needs a display that can handle variable refresh rate input, e.g.,
% a "FreeSync", DisplayPort "Adaptive Sync" or HDMI "VRR" compatible monitor, in
% some cases a "G-Sync" compatible monitor. Also a suitable graphics card with such
% support.
%
% As of December 2018 this function is *highly experimental* and may not work at all,
% or not reliable at all even on nominally suitable system hardware + software
% configurations. Don't rely on it working on any given setup, this is a sneak
% preview of work in progress! Behaviour of the functionality may change at any
% time in backwards incompatible ways without warning!

% This is the numeric constant for kPsychUseFineGrainedOnset:
rc = 2^31;

return;
