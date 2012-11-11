function KbQueueCreate(deviceNumber, keyList)
% KbQueueCreate([deviceNumber, keyList])
%
% The routines KbQueueCreate, KbQueueStart, KbQueueStop, KbQueueCheck
%  KbQueueWait, KbQueueFlush and KbQueueRelease provide replacments for
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
%     2) If many keys are pressed very quickly in succession, it is at least
%        theoretically possible for the queue to fill more quickly than
%        it can be emptied, losing key events temporarily while filled to 
%        capacity. The queue holds up to thirty events, and events are
%        constantly being removed from the queue and processed, so this is
%        unlikely to be a problem in actual use.
%
% Only a single device can be monitored at any given time. The deviceNumber can 
%  be specified only in the call to KbQueueCreate. The other routines then 
%  relate to that specified device. If deviceNumber is not specified, the first 
%  device is the default (like KbCheck). If KbQueueCreate has not been called 
%  first, the other routines will generate an error message. Likewise, if 
%  KbQueueRelease has been called more recently than KbQueueCreate, the other 
%  routines will generate error messages.
%
% It is acceptable to call KbQueueCreate at any time (e.g., to switch to a new
%  device or to change the list of queued keys) without calling KbQueueRelease.
%
%  KbQueueCreate([deviceNumber][, keyList])
%      Creates the queue for the specified (or default) device number
%      If the device number is less than zero, the default device is used.
%      keyList is an optional 256-length vector of doubles (not logicals)
%      with each element corresponding to a particular key (use KbName
%      to map between keys and their positions). If the double value
%      corresponding to a particular key is zero, events for that key
%      are not added to the queue and will not be reported.
%      No events are delivered to the queue until KbQueueStart or 
%      KbQueueWait is called.
%      KbQueueCreate can be called again at any time. The function can also
%      treat other HID devices with buttons or keys as if they are
%      keyboards. E.g., it can also record button state of a mouse, a
%      joystick or a gamepad.
%			
%  KbQueueStart()
%      Starts delivering keyboard events from the specified device to the 
%        queue.
%			
%  KbQueueStop()
%      Stops delivery of new keyboard events from the specified device to 
%        the queue.
%      Data regarding events already queued is not cleared and can be 
%        recovered by KbQueueCheck
%
% [pressed, firstPress, firstRelease, lastPress, lastRelease]=
%   KbQueueCheck()
%      Obtains data about keypresses on the specified device since the 
%        most recent call to this routine, KbQueueStart, KbQueuWait
%      Clears all scored events, but unscored events that are still being
%        processsed may remain in the queue
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
%		ptb_kbcheck_disabledKeys (see "help DisableKeysForKbCheck"), will
%       not caused pressed to return as true and will be zeroed out in the
%       returned arrays. However, a better alternative is to specify a
%       keyList arguement to KbQueueCreate. 
%
% secs=KbQueueWait()
%      Waits for any key to be pressed and returns the time of the press
%
%      KbQueueFlush should be called immediately prior to this function
%      (unless the queue has just been created and started) to clear any 
%      prior events.
%
%      Note that this command will not respond to any keys that were 
%       inactivated by using the keyList argument to KbQueueCreate.
%
%      Since KbQueueWait is implemented as a looping call to
%       KbQueueCheck, it will not respond to any key codes stored in
%       the global variable ptb_kbcheck_disabledKeys
%       (see "help DisableKeysForKbCheck")
%
% KbQueueFlush()
%      Removes all unprocessed events from the queue and zeros out any
%       already scored events.
%
% KbQueueRelease()
%      Releases queue-associated resources; once called, KbQueueCreate
%        must be invoked before using any of the other routines
%
%      This routine is called automatically at clean-up (e.g., when 
%        'clear mex' is invoked and can be omitted expense of keeping 
%        memory allocated and an additional thread running unnecesarily
%
% Note that any keyboard typing used to invoke KbQueue commands will be
%  recorded. This would include the release of the carriage return used
%  to execute KbQueueStart and the keys pressed and released to invoke 
%  KbQueueCheck
% _________________________________________________________________________
%
% See also: KbQueueCreate, KbQueueStart, KbQueueStop, KbQueueCheck,
%            KbQueueWait, KbQueueFlush, KbQueueRelease

% 8/19/07    rpw  Wrote it.
% 8/23/07    rpw  Modifications to add KbQueueFlush

% Requires Mac OS X 10.3 or later. We sort this out on the first call 
% and then store the result in macosrecent for subsequent calls

persistent macosxrecent;

if nargin < 1
    deviceNumber = [];
end

% Try to reserve keyboard queue for 'deviceNumber' for our exclusive use:
if ~KbQueueReserve(1, 2, deviceNumber)
    error('Keyboard queue for device %i already in use by GetChar() et al. Use of GetChar and keyboard queues is mutually exclusive!', deviceNumber);
end

if isempty(macosxrecent)
    macosxrecent = IsOSX;
    LoadPsychHID;
end

if nargin == 2
  PsychHID('KbQueueCreate', deviceNumber, keyList);
elseif nargin == 1
  PsychHID('KbQueueCreate', deviceNumber);
elseif nargin == 0
  PsychHID('KbQueueCreate');
elseif nargin > 2
  error('Too many arguments supplied to KbQueueCreate'); 
end
