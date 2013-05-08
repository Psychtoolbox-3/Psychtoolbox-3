function err=DaqReset(OldDaqIndex)
% err=DaqReset(DeviceIndex)
% Assuming that something's wrong with the USB-1208FS or our communication
% with it, we re-enumerate in order to re-establish communication. Then we
% send the reset command to ask the USB-1208FS to reset itself. Then we
% re-enumerate again to re-establish communication once more.
% 
% To avoid problems caused by CLEAR PsychHID, we recommend that (if you're using
% a 1208FS), instead of calling DaqReset, you unplug and reinsert the USB cable
% of your USB-1208FS and quit and restart MATLAB. In Denis' experience that
% combination always restores normal communication.  If you are using a 1608FS,
% keep reading...
% 
% This function calls "clear PsychHID" twice, and yet I still frequently found
% that I needed to run that command again in order for communication to be
% properly established again.  With a USB-1608FS, Matlab 2007b, and Leopard, I
% found that I didn't have the problems Denis seemed to have.  But what I did
% have was a problem with PsychHID not finding all of the interfaces when
% devices were enumerated.  Running this function (followed by an additional
% "clear PsychHID" command) worked for me, so my recommendation for that case is
% the opposite of Denis'.  I never needed to re-start Matlab or unplug the
% device to get my problems solved.  So I recommend you just run this command,
% then run "clear PsychHID", then try "daq=DaqFind" or "daqs=DaqDeviceIndex"
% (the latter if you have more than one A/D converter built by Measurement
% Computing.  -- mpr 
%
% On Snow Leopard. Matlab R2010a, I found I could re-establish communication
% with an unresponsive 1208FS with the calls above (that is, calling DaqReset
% followed by "clear PsychHID") -- sdv
% 
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest, DaqFind, 
% DaqDeviceIndex.

% 4/15/05 dgp Wrote it.
% 1/9/08  mpr tweaked it for use with 1608FS
% 6/30/2010 sdv fixed error if any USB devices were shorter than 10 characters

fprintf('Resetting USB-1x08FS.\n');
clear PsychHID; % flush current enumeration  (list of devices)
clear PsychHIDDAQS;
devices=PsychHIDDAQS; % enumerate again
daq=[];
for k=1:length(devices)
  if length(devices(k).product)>=10
   if all(devices(k).product([4:5 7:10]) == '-108FS')
    if isempty(daq)
      daq=k;
    elseif ~streq(devices(k).serialNumber,devices(daq(end)).serialNumber)
      daq(end+1)=k;
    end %isempty
   end %all devices
  end %length
end
if isempty(daq)
  error('Sorry, couldn''t find a USB-1x08FS.');
end
if ~any(ismember(OldDaqIndex,daq))
  NewDaqIndex=daq(1);
else
  NewDaqIndex=OldDaqIndex;
end
% Reset. Ask the USB-1208FS to reset its USB interface.
% CAUTION: Immediately after RESET, all commands fail, returning error
% messages saying the command is unsupported (0xE00002C7) or the device is
% not responding (0xE00002ED) or not attached (0xE00002D9). To restore
% communication we must flush the current enumeration and re-enumerate the
% HID-compliant devices.
% fprintf('Reset. Ask the USB-1208FS to reset its USB interface.\n');
err=Reset(NewDaqIndex);
% fprintf('(Reestablishing communication: Flushing current enumeration. ');
clear PsychHID; % flush current enumeration  (list of devices)
clear PsychHIDDAQS;

% fprintf('Re-enumerating. ... ');
devices=PsychHIDDAQS; % enumerate again

% Not sure what the point of this next bit is since we don't do anything with
% the information we acquire, but it seems harmless enough so I leave it in. --
% mpr
daq=[];
for k=1:length(devices)
 if length(devices(k).product)>=10,
  if all(devices(k).product([4:5 7:10]) == '-108FS')
    if isempty(daq)
      daq=k;
    elseif ~streq(devices(k).serialNumber,devices(daq(end)).serialNumber)
      daq(end+1)=k;
    end
  end
 end
end
% fprintf('Done.)\n');
return

function err=Reset(DaqIndex)
% err=Reset(DeviceIndex)
% USB-1208FS: Ask device to reset its USB interface.
err=PsychHID('SetReport',DaqIndex,2,65,uint8(65)); % Reset
if err.n
  fprintf('Reset SetReport error 0x%s. %s: %s\n',hexstr(err.n),err.name,err.description);
end
return

