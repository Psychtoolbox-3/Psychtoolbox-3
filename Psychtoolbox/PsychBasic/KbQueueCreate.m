function KbQueueCreate(deviceNumber, keyList, numValuators, numSlots, flags, win)
% KbQueueCreate([deviceNumber][, keyList][, numValuators=0][, numSlots=10000][, flags=0][, windowHandle=0])
%
% The routines KbQueueCreate, KbQueueStart, KbQueueStop, KbQueueCheck
%  KbQueueWait, KbQueueFlush and KbQueueRelease provide replacements for
%  KbCheck and KbWait, providing the following advantages:
%
%     1) Brief key presses that would be missed by KbCheck or KbWait
%        are reliably detected
%     2) The times of key presses are recorded more accurately
%     3) The times of key releases are also recorded
%
% Limitations:
%
%     1) If a key is pressed multiple times before KbQueueCheck is called,
%        only the times of the first and last presses and releases of that
%        key can be recovered (this has no effect on other keys)
%     2) If many keys are pressed very quickly in succession on OSX, it is
%        at least theoretically possible for the queue to fill more quickly
%        than it can be emptied, losing key events temporarily while filled
%        to capacity. The queue holds up to thirty events, and events are
%        constantly being removed from the queue and processed, so this is
%        unlikely to be a problem in actual use.
%
%  The deviceNumber of a HID input device can be specified as 'deviceNumber'.
%  Allowable devices are "keyboard like" devices, e.g., Keyboards, Keypads,
%  Mouse, Joystick, Gamepad, Touchpad, etc. - Stuff that has buttons or keys
%  to press. This way a mouse or joysticks buttons can be used as response
%  devices, ignoring mouse movements etc. If deviceNumber is not specified, the first 
%  device is the default (like KbCheck). If KbQueueCreate has not been called 
%  first, the other routines will generate an error message. Likewise, if 
%  KbQueueRelease has been called more recently than KbQueueCreate, the other 
%  routines will generate error messages.
%
% It is acceptable to call KbQueueCreate at any time (e.g., to switch to a new
%  device or to change the list of queued keys) without calling KbQueueRelease.
%
%  KbQueueCreate([deviceNumber][, keyList][, numValuators=0][, numSlots=10000][, flags=0][, windowHandle=0])
%      Creates the queue for the specified (or default) device number
%      If the device number is less than zero, the default device is used.
%
%      'keyList' is an optional 256-length vector of doubles (not logicals)
%      with each element corresponding to a particular key (use KbName
%      to map between keys and their positions). If the double value
%      corresponding to a particular key is zero, events for that key
%      are not added to the queue and will not be reported.
%
%      'numValuators' is an optional maximum number of additional values to report.
%      It defaults to zero. For values greater than zero, if the selected type
%      of input device supports this, and if the operating system supports this,
%      additional info will be recorded. For pointing devices like mice, the mouse
%      position may be reported, for joysticks the state of their various axis, etc.
%      For mouse/joystick/pointing device position reporting numValuators must be
%      at least 2. On touch devices, a value of 2 treats them like a mouse, a value
%      >= 4 treats them as touchscreen. However, for touch devices there are separate
%      functions like TouchQueueCreate, TouchEventGet etc. for convenient handling.
%      See "help KbEventGet" for how to retrieve potential additionally recorded
%      information.
%
%      'numSlots' defines how many events the event buffer can store. If a script
%      does not periodically remove events via KbEventGet() or KbEventFlush(), the
%      buffer will fill up, and once 'numSlots' elements are stored, it will stop
%      recording new events. 10000 elements capacity is the default, which may be
%      too little if you use 'numValuators' > 0 to store dynamic (motion) data like
%      mouse movements or touchscreen input, which can be generated at rates of
%      multiple hundred events per second of data collection.
%
%      'flags' defines special modes of operation for the queue. These are OS
%      specific, see "PsychHID KbQueueCreate?" for an up to date list of supported
%      flags. In general, you don't need these.
%
%      'windowHandle' defines the optional onscreen window handle of an associated
%      onscreen window. By default, input is taken for all windows and screens.
%      This argument is silently ignored on systems other than Linux/X11 at the moment.
%
%      No events are delivered to the queue until KbQueueStart or 
%      KbQueueWait is called.
%      KbQueueCreate can be called again at any time. The function can also
%      treat other HID devices with buttons or keys as if they are
%      keyboards. E.g., it can also record button state of a mouse, a
%      joystick or a gamepad.
%
%  KbQueueStart([deviceNumber])
%      Starts delivering keyboard events from the specified device to the 
%      queue.
%
%  KbQueueStop([deviceNumber])
%      Stops delivery of new keyboard events from the specified device to 
%      the queue.
%      Data regarding events already queued is not cleared and can be 
%      recovered by KbQueueCheck
%
% [pressed, firstPress, firstRelease, lastPress, lastRelease]=
%   KbQueueCheck([deviceNumber])
%      Obtains data about keypresses on the specified device since the 
%      most recent call to this routine, KbQueueStart, KbQueueWait
%      Clears all scored events, but unscored events that are still being
%      processsed may remain in the queue
%
%      pressed: a boolean indicating whether a key has been pressed
%
%      firstPress: an array indicating the time that each key was first
%        pressed since the most recent call to KbQueueCheck or KbQueueStart
%
%      firstRelease: an array indicating the time that each key was first
%        released since the most recent call to KbQueueCheck or KbQueueStart
%
%      lastPress: an array indicating the most recent time that each key was
%        pressed since the most recent call to KbQueueCheck or KbQueueStart
%
%      lastRelease: an array indicating the most recent time that each key
%         was released since the most recent call to KbQueueCheck or 
%         KbQueueStart
%
%     For firstPress, firstRelease, lastPress and lastRelease, a time value
%       of zero indicates that no event for the corresponding key was
%       detected since the most recent call to KbQueueCheck or KbQueueStart
%
%     To identify specific keys, use KbName (e.g., KbName(firstPress)) to
%       generate a list of the keys for which the events occurred
%
%     For compatibility with KbCheck, any key codes stored in
%     ptb_kbcheck_disabledKeys (see "help DisableKeysForKbCheck"), will
%     not cause pressed to return as true and will be zeroed out in the
%     returned arrays. However, a better alternative is to specify a
%     keyList arguement to KbQueueCreate. 
%
% secs=KbQueueWait([deviceNumber])
%      Waits for any key to be pressed and returns the time of the press.
%
%      KbQueueFlush should be called immediately prior to this function
%      (unless the queue has just been created and started) to clear any 
%      prior events.
%
%      Note that this command will not respond to any keys that were 
%      inactivated by using the keyList argument to KbQueueCreate.
%
%      Since KbQueueWait is implemented as a looping call to
%      KbQueueCheck, it will not respond to any key codes stored in
%      the global variable ptb_kbcheck_disabledKeys
%      (see "help DisableKeysForKbCheck")
%
% KbQueueFlush([deviceNumber])
%      Removes all unprocessed events from the queue and zeros out any
%      already scored events.
%
% KbQueueRelease([deviceNumber])
%      Releases queue-associated resources; once called, KbQueueCreate
%      must be invoked before using any of the other routines
%
%      This routine is called automatically at clean-up (e.g., when 
%      'clear mex' is invoked and can be omitted expense of keeping 
%      memory allocated and an additional thread running unnecesarily
%
% Note that any keyboard typing used to invoke KbQueue commands will be
% recorded. This would include the release of the carriage return used
% to execute KbQueueStart and the keys pressed and released to invoke 
% KbQueueCheck
% _________________________________________________________________________
%
% See also: KbQueueCreate, KbQueueStart, KbQueueStop, KbQueueCheck,
%           KbQueueWait, KbQueueFlush, KbQueueRelease

% 8/19/07    rpw  Wrote it.
% 8/23/07    rpw  Modifications to add KbQueueFlush

persistent macosxrecent;

if nargin < 1
  deviceNumber = [];
end

% Try to reserve keyboard queue for 'deviceNumber' for our exclusive use:
if ~KbQueueReserve(1, 2, deviceNumber)
  if isempty(deviceNumber)
      deviceNumber = NaN;
  end
  error('Keyboard queue for device %i already in use by GetChar() et al. Use of GetChar and keyboard queues is mutually exclusive!', deviceNumber);
end

if isempty(macosxrecent)
  macosxrecent = IsOSX;
  LoadPsychHID;
end

if nargin >= 6
  if ~isempty(win)
    % Onscreen window handle given. Validate and map to windowing system handle:
    if Screen('WindowKind', win) == 1 
      % Onscreen window handle of open onscreen window: Map to X-Window:
      winfo = Screen('GetWindowInfo', win);
      win = winfo.SysWindowHandle;
    elseif ismember(win, Screen('Screens'))
      % screenid given: Pass in as X-Screen id + 1:
      win = win + 1;
    end
  end
end

if nargin == 6
  PsychHID('KbQueueCreate', deviceNumber, keyList, numValuators, numSlots, flags, win);
elseif nargin == 5
  PsychHID('KbQueueCreate', deviceNumber, keyList, numValuators, numSlots, flags);
elseif nargin == 4
  PsychHID('KbQueueCreate', deviceNumber, keyList, numValuators, numSlots);
elseif nargin == 3
  PsychHID('KbQueueCreate', deviceNumber, keyList, numValuators);
elseif nargin == 2
  PsychHID('KbQueueCreate', deviceNumber, keyList);
elseif nargin == 1
  PsychHID('KbQueueCreate', deviceNumber);
elseif nargin == 0
  PsychHID('KbQueueCreate');
elseif nargin > 4
  error('Too many arguments supplied to KbQueueCreate'); 
end
