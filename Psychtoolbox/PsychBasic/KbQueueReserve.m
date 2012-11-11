function isReserved = KbQueueReserve(action, actor, deviceIndex)
% Reserve the keyboard queue of the default keyboard for use by the
% alternative GetChar() implementation. This is for internal use only!
%
% 'action': 1 = Try to reserve, 2 = Try to release, 3 = Query reservation.
% 'actor' : For whom to reserve/release/query ownership: 1 = GetChar, 2 = Usercode
%
% The function will reserve or release the queue on behalf of 'actor' if it
% isn't already reserved for another actor.
%
% The function returns 1 if the queue is now reserved to 'actor', 0
% otherwise.
%

% History:
% 23.10.2012  mk  Written.

% Store for whom the default queue is reserved:
persistent reservedFor;
persistent defaultKbDevice;

if isempty(reservedFor)
    % Initially not reserved for anybody:
    reservedFor = 0;

    % Get deviceIndex of default keyboard device for KbQueues:
    LoadPsychHID;
    defaultKbDevice = PsychHID('Devices', -1);
end

% Special case handling for Linux and Windows:
if ~IsOSX && ~isempty(deviceIndex) && (deviceIndex ~= defaultKbDevice)
    % On non-OSX, all non-default-keyboard queues are always reserved for
    % usercode, as only the default keyboard queue (aka empty deviceIndex)
    % matters for GetChar:
    if actor == 2
        isReserved = 1;
    else
        isReserved = 0;
    end
    
    return;
end

% Either OSX, or MS-Windows/Linux on default keyboard device. There's only
% one such deviceIndex zero/default on non-OSX, and only one keyboard queue in
% total on OSX, irrespective of deviceIndex. Therefore a simple variable is enough
% to keep reservation status for that one queue.

% Reserve request?
if action == 1
    % If it is already reserved for us, or not reserved to anybody, then we
    % can reserve it for us:
    if (reservedFor == 0) || (reservedFor == actor)
        
        if (actor == 1) && IsOSX && ~Is64Bit
            % Special case: Tried to use KbQueue for GetChar() et al. on
            % 32-Bit OSX runtime, which is unsupported due to limitations
            % in PsychHID:
            error('GetChar/CharAvail/FlushEvents/ListenChar are not supported on 32-Bit Matlab in -nojvm mode on OSX.');
        end
        
        reservedFor = actor;
    end
end

% Release request?
if action == 2
    % If it is reserved for us, or not reserved to anybody, then we
    % can safely release it, so it does not belong to anybody:
    if (reservedFor == 0) || (reservedFor == actor)
        reservedFor = 0;
    end
end

% Return 1 = True, if queue is now reserved for us, 0 = False otherwise.
if reservedFor == actor
    % Reserved for us:
    isReserved = 1;
else
    % Not available for us:
    isReserved = 0;
end

return;
