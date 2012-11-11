function secs = KbTriggerWait(keyCode, deviceNumber)
% secs = KbTriggerWait(keyCode, [deviceNumber])
%
% Waits until one or more trigger keys have been pressed and returns the
% time of the first key press in seconds. The keyCode argument can be a
% vector of key indices. For example, to check for the 't' key as the
% trigger use KbTriggerWait(KbName('t')), to check for 't' and the escape
% key, use KbTriggerWait([KbName('t'), KbName('ESCAPE')]).
%
% You cannot use KbTriggerWait while a queue created by KbQueueCreate
% exists. To shut down such a queue, use KbQueueRelease.
%
% On Matlab versions older than R2007a on MS-Windows, this function simply
% serves as a convenient substitute for using KbCheck to detect the
% trigger of interest.
%
% This function should allow triggers to be reliably detected from devices
% that only briefly report that the key is down. KbCheck is not reliable
% with such devices because it may not poll often enough to detect the
% key down state.
%
% KbTriggerWait uses the PsychHID function, a general purpose function for
% reading from the Human Interface Device (HID) class of USB devices.
% Unlike KbCheck, it starts a queue that receives keyboard events
% regarding the trigger key (and no other keys) and then polls this queue
% (rather than the current key status) periodically. In theory, this
% should also provide more accurate reporting of the time of the
% triggering keypress. However, if multiple trigger events have occurred
% since last polled, it is possible that the timestamp of the earliest
% of these will have already rotated out of the limited capacity (eight
% events) queue. In this case, the time of the earliest event remaining
% in the queue is reported. Since the polling frequency is the same as
% KbCheck, it should be more accurate on average with regard to timing,
% even when the timestamps of the earliest events have been lost due to
% queue overflow.
%
% KbTriggerWait tests the first USB-HID keyboard device by default.
% Optionally, you can pass in a 'deviceNumber' to test a different keyboard
% if multiple keyboards are connected to your machine. The function also
% allows to wait for button presses on keypads, mice or other HID devices
% with buttons or keys.
%
% Passing a deviceNumber of -1 will NOT cause all keyboards to be detected
%
% One disadvantage of this function is that it renders Matlab relatively
% unresponsive to Ctrl-C interrupts. KbQueueWait is a better option in
% this regard, but more complicated to use.
% _________________________________________________________________________
%
% See also: KbQueueWait KbCheck, KbWait, GetChar, CharAvail, KbDemo.

% 8/10/07    rpw  Wrote it.
% 8/21/07	 rpw  Added comments about KbQueueWait as alternative
% 8/23/07    rpw  Added warning about incompatibility with KbQueueCreate, et al.
% 5/14/12    mk   Add new OSX path, MS-Windows pre R2007a path, other
%                 tweaks.

persistent macosxrecent;
if isempty(macosxrecent)
    macosxrecent = IsOSX;
    LoadPsychHID;
end

% OSX? We no longer use PsychHID('KbTriggerWait'). Instead we emulate via
% the KbQueueXXX interface. This has two advantages:
% a) It works on 64-Bit OSX despite the incompatibilities and additional
%    brain-damage added to the low-level interface by our friendly iPhone
%    company.
% b) It allows to wait for multiple keys, just as on Linux and Windows.
%
% Timestamp precision won't suffer, as past improvements to PsychHID's
% KbQueue routines now query timestamps from the OS itself for high
% precision.
%
if macosxrecent
    if nargin==2
        % Nothing to do.
    elseif nargin == 1
        deviceNumber = [];
    elseif nargin == 0
        error('Trigger key code(s) must be specified in KbTriggerWait');
    elseif nargin > 2
        error('Too many arguments supplied to KbTriggerWait');
    end
    
    % Emulate KbTriggerWait via KbQueueWait:
    keyCodes = zeros(1, 256);
    keyCodes(keyCode) = 1;
    KbQueueCreate(deviceNumber, keyCodes);
    secs = KbQueueWait(deviceNumber);
    KbQueueStop(deviceNumber);
    KbQueueRelease(deviceNumber);
    return;
end

if nargin < 2
    deviceNumber = [];
end

% Try to reserve keyboard queue for 'deviceNumber' for our exclusive use:
if ~KbQueueReserve(1, 2, deviceNumber)
    error('Keyboard queue for device %i already in use by GetChar() et al. Use of GetChar and keyboard queues is mutually exclusive!', deviceNumber);
end

if nargin==2
    [secs]= PsychHID('KbTriggerWait', keyCode, deviceNumber);
elseif nargin == 1
    [secs]= PsychHID('KbTriggerWait', keyCode);
elseif nargin == 0
    error('Trigger key code must be specified in KbTriggerWait');
elseif nargin > 2
    error('Too many arguments supplied to KbTriggerWait');
end

% Try to release keyboard queue for 'deviceIndex' from our exclusive use:
KbQueueReserve(2, 2, deviceNumber);
