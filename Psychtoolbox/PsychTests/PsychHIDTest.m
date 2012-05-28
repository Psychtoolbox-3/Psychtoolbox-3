function PsychHIDTest
% PsychHIDTest
%
% PsychHIDTest exercises the PsychHID mex file. We list all the HID
% devices. We read from the keyboard and mouse. We flicker the keyboard
% LEDs.
% 
% NOT RESPONDING? If PsychHID is not responding, e.g. after unplugging and
% re-plugging the USB connector, try quitting and restarting MATLAB. We
% find that this reliably restores normal communication. 
% 
% web http://psychtoolbox.org/usb.html -browser;
% See also DaqTest, PsychHID, PsychHardware.

% 3/19/05 dgp Denis Pelli wrote it.
% 4/18/05 dgp Flush stale reports.
% 8/26/05 dgp Incorporated bug fixes (for compatibility with Mac OS X Tiger) 
%             suggested by Jochen Laubrock <laubrock@rz.uni-potsdam.de>
%             and Maria Mckinley <parody@u.washington.edu>. 
%             The reported number of outputs of the USB-1208FS has changed in Tiger.
%             http://groups.yahoo.com/group/psychtoolbox/message/3610
%             http://groups.yahoo.com/group/psychtoolbox/message/3614
% 5/14/12  mk Cleanup and improve.

fprintf('TestPsychHID\n');
fprintf('Making a list of all your HID-compliant devices. ...');
devices=PsychHID('Devices');
fprintf('\n\nYou have %d HID-compliant devices:\n',length(devices));
for di=1:length(devices)
    d=devices(di);
    s=sprintf('device %d: %s, %s, %s',di,d.usageName,d.manufacturer,d.product);
    s=sprintf('%s, %d inputs, %d outputs',s,d.inputs,d.outputs);
    if ~isempty(d.serialNumber)
        s=sprintf('%s, serialNumber %s',s,d.serialNumber);
    end
    fprintf('%s\n',s);
end
fprintf('\n');

% USB-1208FS
% daq=[];
% for di=1:length(devices)
%     if (streq(devices(di).product,'PMD-1208FS') || streq(devices(di).product,'USB-1208FS')) && devices(di).outputs>=70
%         % Select the main interface to represent the whole unit
%         daq(end+1)=di; %#ok<AGROW>
%     end
% end
daq = DaqDeviceIndex;

switch length(daq)
    case 0;
    case 1,
        fprintf('You have a USB-1208FS. It appears as four devices in the table above.\n');
    otherwise,
        fprintf('You have %d USB-1208FS boxes. Each appears as four devices in the table above.\n',length(daq));
end

% Keyboard.
fprintf('\n');
keyboard = GetKeyboardIndices;
switch length(keyboard)
    case 0;
    case 1,
        fprintf('Test keyboard.\n');
    otherwise,
        fprintf('Test %d keyboards.\n',length(keyboard));
end
for di=keyboard
    fprintf('device %d: %s, %s\n',di,devices(di).manufacturer,devices(di).product);
    if devices(di).outputs>0
        fprintf('Now flickering the LEDs on your keyboard. ... ');
        for i=1:16*6 % for 6 s
            WaitSecs(1/16);
            err=PsychHID('SetReport',di,2,0,uint8(255*mod(i,2)));
            if err.n
                break
            end
        end
        if err.n
            fprintf('\nPsychHID: SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        else
            fprintf('Success!!\n');
        end
    end
    r=1;
    while ~isempty(r)
        % flush any old reports
        [r,err]=PsychHID('GetReport',di,1,0,8);
    end
    fprintf('Now reading your keyboard. Press any key to continue. ... ');
    WaitSecs(0.1); % Wait a moment for key release.
    keydata=[];
    timeout=GetSecs+5;
    while ~any(keydata) && GetSecs<timeout
        [r,err]=PsychHID('GetReport',di,1,0,8);
        if err.n
            fprintf('\nPsychHID: GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        end
        if length(r)>2
            keydata=r([1 3:end]); % Omit second byte, which has no data.
        end
    end
    if any(keydata)
        fprintf('Thanks.\n');
    else
        fprintf('\nI gave up waiting after 5 s.\n');
    end
    PsychHID('ReceiveReportsStop',di);
end

% Mouse.
fprintf('\n');
mouse = GetMouseIndices;
switch length(mouse)
    case 0;
    case 1,
        fprintf('Test mouse.\n');
    otherwise,
        fprintf('Test %d mice.\n',length(mouse));
end
for di=mouse
    r=1;
    while ~isempty(r)
        % flush any old reports
        r = PsychHID('GetReport',di,1,0,8);
    end
    fprintf('device %d: %s, %s\n',di,devices(di).manufacturer,devices(di).product);
    fprintf('Now reading your mouse. Move the mouse to continue. ... ');
    mousedata=[];
    timeout=GetSecs+5;
    while GetSecs<timeout && ~any(mousedata)
        [r,err]=PsychHID('GetReport',di,1,0,3);
        if err.n<0
            fprintf('\nPsychHID: GetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
        end
        if ~isempty(r)
            mousedata=r(2:3);
        end
    end
    if any(mousedata)
        fprintf('Thanks.\n');
    else
        fprintf('\nI gave up waiting after 5 s.\n');
    end
    PsychHID('ReceiveReportsStop',di);
end

return
% The code below prints out lots of information, but I
% haven't found it to be particularly interesting. Comment out the
% return statement, immediately above, to get the extra info.

% elements
fprintf('\n');
for di=1:length(devices)
    fprintf('device %d ******************\n',di);
    devices(di)
    numberOfElements=PsychHID('NumElements',di);
    fprintf('numberOfElements %d\n',numberOfElements);
    elements=PsychHID('Elements',di);
    elements(1)
    for e=1:length(elements)
        %         state=PsychHID('RawState',di,e);
        %         calibratedState=PsychHID('CalibratedState',di,e);
        % fprintf('%d,elementState raw %d calibrated %d\n',e,state,calibratedState);
    end
    numberOfCollections=PsychHID('NumCollections',di);
    fprintf('numberOfCollections %d\n',numberOfCollections);
    collections=PsychHID('Collections',di);
    collections(1)
end

function str=hexstr(n)
% str=hexstr(n)
% Convert any number to a hex string representing the lower 32 bits,
% emulating the C format %lx for a long. This works with both positive and
% negative numbers, unlike the MATLAB format %x (and %tx and %bx) which
% can't deal with negative numbers.
%
% fprintf('Error 0x%s.\n',hexstr(err));
%
% See also DEC2HEX, FPRINTF.

% 4/24/05 dgp Wrote it.

h=dec2hex(n+2^33);
str=h(end-7:end);
