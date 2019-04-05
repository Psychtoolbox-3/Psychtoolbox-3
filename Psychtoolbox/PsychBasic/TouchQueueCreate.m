function TouchQueueCreate(win, deviceNumber, numSlots, numValuators, keyList, flags)
% TouchQueueCreate(windowHandle, deviceNumber [, numSlots=100000][, numValuators=auto][, keyList=all][, flags=0])
%
% Create a touch queue for receiving touch input from touch input devices
% like touchscreens, tablets, touch surfaces, or touchpads.
%
% 'windowHandle' Handle for a Screen() onscreen window for which touch
% should be received.
%
% 'deviceNumber' Handle for the touch input device. GetTouchDeviceIndices()
% allows to enumerate touch devices.
%
% 'numSlots' Number of input slots for storing touch events. You must empty
% the queue frequently enough, so no more than 'numSlots' events collect in
% the queue, or the queue will stop storing new events and drop information.
%
% 'numValuators' Number of optional touch device properties to store per
% event. Defaults to what the device can provide.
%
% 'keyList' If the touch queue is also used to receive button input from
% physical buttons (or maybe virtual buttons?) from the device, this provides
% the list of buttons to accept. See KbQeueCreate for explanation. By default
% all buttons are accepted.
%
% 'flags' Optional flags. This defaults to 0 and is so far unused for touch input.
%
% Once a queue is created its touch data collection can be started via
% TouchQueueStart(), stopped via TouchQueueStop(), cleared via TouchQueueFlush(),
% asked for current number of pending events via TouchEventAvail() and events
% can be fetched via TouchEventGet().
%

% History:
% 1-Oct-2017  mk  Written.

if ~IsLinux
  error('This function is currently only supported on Linux/X11');
end

if nargin < 1 || isempty(win) || Screen('WindowKind', win) ~= 1
  error('Required windowHandle missing.');
end

if nargin < 2 || isempty(deviceNumber)
  error('Required deviceNumber missing.');
end

if nargin < 3 || isempty(numSlots)
  numSlots = 100000;
else
  if numSlots < 4
    error('numSlots must be at least 4.');
  end
end

if nargin < 4
  numValuators = [];
else
  if numValuators < 4
    error('numValuators must be at least 4.');
  end
end

if nargin < 5
  keyList = [];
else
  if ~ismember(length(keyList), [0, 256])
    error('keyList must be either empty, or have exactly 256 entries.');
  end
end

if nargin < 6
  flags = [];
end

[touchIndices, productNames, allInfo] = GetTouchDeviceIndices;
if ~ismember(deviceNumber, touchIndices)
  error('deviceNumber does not refer to a touch input device.');
end

% Get dev info for touch device:
idx = find(touchIndices == deviceNumber);
dev = allInfo{idx};

if isempty(numValuators)
  if IsLinux
    % On Linux/X11 this is stored in axes, but we need at least numValuators of
    % 4, as otherwise touch input support won't get enabled by PsychHID():
    numValuators = max(dev.axes, 4);
  else
    % On MS-Windows, 5 would be a good guess:
    numValuators = 5;
  end
end

% Ok, and now we just call KbQeueCreate, which we wrap here atm.:
KbQueueCreate(deviceNumber, keyList, numValuators, numSlots, flags, win);

return;
