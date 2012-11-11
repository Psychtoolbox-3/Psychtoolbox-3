function nflushed = KbEventFlush(deviceIndex)
% nflushed = KbEventFlush([deviceIndex])
%
% Flush event buffer of a keyboard queue. This removes all stored events
% from the keyboard event buffer of a given keyboard queue. It returns the
% number of removed events in the optional return argument 'nflushed'.
%
% By default, the event buffer of the default keyboard queue is emptied,
% but you can specify 'deviceIndex' to flush the buffer of the queue
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

nflushed = PsychHID('KbQueueFlush', deviceIndex, 2);

return;
