function [event, nremaining] = KbEventGet(deviceIndex, maxWaitTimeSecs)
% [event, nremaining] = KbEventGet([deviceIndex][, maxWaitTimeSecs=0])
%
% Return oldest pending event, if any, in return argument 'event', and the
% remaining number of recorded events in the event buffer of a keyboard
% queue in the return argument 'nremaining'.
%
% By default, the event buffer of the default keyboard queue is checked,
% but you can specify 'deviceIndex' to check the buffer of the queue
% associated with 'deviceIndex'.
%
% KbEventGet() will wait up to 'maxWaitTimeSecs' seconds for at least one
% event to show up before it gives up. By default, it doesn't wait but just
% gives up if there aren't any events queued at time of invocation.
%
% 'event' is either empty if there aren't any events available, or it is a
% struct with information about the keyboard event. The returned event
% struct currently contains the following fields:
%
% 'Keycode' = The KbCheck / KbName style keycode of the key or button that
%             triggered this event.
%
% 'Time' = The GetSecs time of when the event was received.
%
% 'Pressed' = 1 for a key press event, 0 for a key release event.
%
% 'CookedKey' = Keycode translated into a GetChar() style ASCII character code.
% Or zero if key does not have a corresponding character. Or -1 if mapping
% is unsupported for given event. This does not yet work correctly on OSX.
%
% Keyboard event buffers are a different way to access the information
% collected by keyboard queues. Before you can use an event buffer you
% always must create its "parent keyboard queue" via KbQueueCreate() and
% call KbQueueStart() to enable key event recording. See "help
% KbQueueCreate" etc. on how to do this.
%
% _________________________________________________________________________
%
% See also: KbQueueCreate, KbQueueStart, KbQueueStop, KbQueueCheck,
%            KbQueueWait, KbQueueFlush, KbQueueRelease

% 21.5.2012  mk  Wrote it.

if nargin < 1
    deviceIndex = [];
end

% Try to check if keyboard queue for 'deviceIndex' is reserved for our exclusive use:
if ~KbQueueReserve(3, 2, deviceIndex) && KbQueueReserve(3, 1, deviceIndex)
    error('Keyboard queue for device %i already in use by GetChar() et al. Use of GetChar and keyboard queues is mutually exclusive!', deviceIndex);
end

if ~KbQueueReserve(3, 2, deviceIndex)
    error('Keyboard queue for device %i not yet created. KbQueueCreate() it first!\n', deviceIndex);
end

if nargin < 2
    maxWaitTimeSecs = [];
end

[event, nremaining] = PsychHID('KbQueueGetEvent', deviceIndex, maxWaitTimeSecs);

return;
