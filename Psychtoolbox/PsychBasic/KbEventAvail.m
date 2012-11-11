function navail = KbEventAvail(deviceIndex)
% navail = KbEventAvail([deviceIndex])
%
% Return current number of recorded events in the event buffer of a
% keyboard queue in the return argument 'navail'.
%
% By default, the event buffer of the default keyboard queue is checked,
% but you can specify 'deviceIndex' to check the buffer of the queue
% associated with 'deviceIndex'.
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

navail = PsychHID('KbQueueFlush', deviceIndex, 0);

return;
