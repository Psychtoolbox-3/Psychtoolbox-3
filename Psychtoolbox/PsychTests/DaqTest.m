function DaqTest
% DaqTest
%
% DaqTest assesses the Daq Toolbox, which provides communication with a
% particular USB data acquisition device (daq): the USB-1208FS made by
% Measurement Computing (see URL below). This daq costs $150 and offers "50
% kHz" input and output 12-bit sampling of analog voltages (8 in, 2 out)
% and 16 digital i/o lines, with signals brought out to screw terminals.
% ("50 kHz" is a theoretical upper limit: as of 18 April 2005 we attain 2
% kHz.) The USB-1208FS is the size of a wallet and is powered through its
% USB cable. The Daq Toolbox gives you complete control of it from within
% Matlab, via the PsychHID extension.
%
% DaqTest assesses all the functions provided by the USB-1208FS firmware.
% (Some risky tests, below, have been turned off, but can be re-enabled by
% the adventurous user.) DaqFunctions gives a brief description of every
% function.
%
% NOT RESPONDING? If PsychHID is not responding, e.g. after unplugging and 
% re-plugging the USB connector of your device, try quitting and restarting
% MATLAB. We find that this reliably restores normal communication. 
%
% web http://www.http://www.measurementcomputing.com/cbicatalog/directory.asp?dept_id=403 -browser;
% web http://psychtoolbox.org/daq.html -browser;
% See also: Daq, DaqFunctions, DaqPins,PsychHIDTest,PsychHID,PsychHardware,
% DaqDeviceIndex, DaqDIn, DaqDOut, DaqAIn, DaqAOut, DaqAInScan, DaqAOutScan

% 3/21/05 dgp Denis Pelli wrote it.
% 3/28/05 dgp Corrected reportID and reportBytes.
% 3/30/05 dgp More error checking.
% 3/31/05 dgp Added the rest of the functions.
% 4/2/05  dgp The input functions now work, for the first time.
% 4/15/05 dgp Polished.
% 4/15/05 dgp Merged arguments in calls to DaqAOutScan and DaqAInScan.
% 4/25/05 dgp Nearly eliminated crash associated with CLEAR MEX.
% 4/26/05 dgp Fixed DaqAInScan to return sensible results.
% 8/26/05 dgp Incorporated bug fix for compatibility with Mac OS X Tiger 
%             suggested by Maria Mckinley <parody@u.washington.edu>. The reported
%             number of outputs of the USB-1208FS has changed in Tiger.
%             http://groups.yahoo.com/group/psychtoolbox/message/3614
% 11/xx/07 - 
%  1/xx/08 mpr  Instituted various fixes to get this to work for a USB-1608FS
%                   in Leopard
% 05/11/09 mk Add if isempty(FinalStatus) in line 444 and later to catch
%             empty returns from DaqGetStatus.

debugging=0; % Used by dgp to track down a crash in PsychHID.

% Do we have a USB-1208FS daq? Each such box will present itself as four
% HID "devices" sharing the same serial number.  For 1608FS, there may be 5, 6,
% or 7 devices.  Don't know why the variation...  Should be 7...  Most commands
% will work even if only five interfaces are found, but the most important
% commands -- the ones that read digitized analog signals -- require all seven
% interfaces unless you're willing to lose some data.  So we'll try to make sure
% all are found before continuing.
daq=DaqDeviceIndex([],0);

if isempty(daq)
  % check again because user might have run PsychHID('Devices') and then plugged
  % in the daq.  This should fix that.
  clear PsychHID;
  daq = DaqDeviceIndex([],0);
end

devices=PsychHID('Devices');
fprintf('\n');
switch length(daq)
  case 0,
    fprintf('Sorry. Couldn''t find any USB-1208FS, -1408FS, or -1608FS boxes connected to your computer.\n');
    return;
  case 1,
    fprintf('Yay. I found a %s daq!  Serial number is: %s\n',devices(daq).product,devices(daq).serialNumber);
    if strcmp(devices(daq).product(5:6),'16')
      eval(['[' sprintf('AllSNs{%d},',1:(length(devices)-1)) 'AllSNs{' int2str(length(devices)) '}] = deal(devices.serialNumber);']);
      NumInterfaces = length(cell2mat(strfind(AllSNs,devices(daq).serialNumber)));
      if NumInterfaces ~= 7
        ConfirmInfo(sprintf('Found %d interfaces where there should be 7.  Run "help DaqReset" for suggestions; then re-run "DaqTest"',NumInterfaces));
        return;
      end
    end
    numTries = 0;
    TheLEDBlinked = 0;
    while ~TheLEDBlinked
      ConfirmInfo('Direct your attention to the LED on the device and then click the "okay" button here:');
      err=DaqBlinkLED(daq);
      if err.n
        fprintf(['\nUh oh... there was a problem when I told your device to blink its LED.\n' ...
                 'I haven''t seen this before, so I don''t know what''s going on.  This is\n' ...
                 'what PsychHID told me about the problem, though:\n\n']);
        err
        fprintf('\nI''ll let you try to fix the problem and then come back to me.  Bailing now!\n\n');
        return;
      end
      TheLEDBlinked = TwoStateQuery('Did the LED just blink for you?');
      if TheLEDBlinked < 0
        ConfirmInfo('Please click yes or no... do not close the window any other way!');
        TheLEDBlinked = TwoStateQuery('Did the LED blink for you earlier?');
        if TheLEDBlinked < 0
          ConfirmInfo('You do not seem to play well with others.  I am quitting on you!');
          return;
        end
      end
      numTries = numTries+1;
      if numTries > 1 && ~TheLEDBlinked
        ConfirmInfo('It should have blinked.  Since it didn''t and there were no obvious errors I will bail and let you work on it.');
        return;
      end
    end
  case 2,
    fprintf('Yay. You have two daqs!');
    numTries = 0;
    TheLEDBlinked = 0;
    while ~TheLEDBlinked
      if strcmp(devices(daq(1)).product(5:6),'16')
        eval(['[' sprintf('AllSNs{%d},',1:(length(devices)-1)) 'AllSNs{' int2str(length(devices)) '}] = deal(devices.serialNumber);']);
        NumInterfaces = length(cell2mat(strfind(AllSNs,devices(daq(1)).serialNumber)));
        if NumInterfaces ~= 7
          ConfirmInfo(sprintf('Found %d interfaces where there should be 7.  Run "help DaqReset" for suggestions; then re-run "DaqTest"',NumInterfaces));
          return;
        end
      end
      fprintf('\nI''m about to try device with serial number: %s\n',devices(daq(1)).serialNumber);
      ConfirmInfo('Look at the LED on the daq you want to test and then click this button:');
      err=DaqBlinkLED(daq(1));
      if err.n
        fprintf(['\nUh oh... there was a problem when I told your first device to blink its LED.\n' ...
                 'This is the error message I received:\n\n']);
        err

        fprintf('\nI will move on to the second one anyways and hope for the best.\n');
      end

      TheLEDBlinked = TwoStateQuery('Did the LED blink on the device you want to test?');
      if TheLEDBlinked < 0
        fprintf('\nYou should have clicked "Yes" or "No".  I am taking your response as a "No".\n');
        TheLEDBlinked = 0;
      end
      if TheLEDBlinked
        daq = daq(1);
        break;
      else
        if strcmp(devices(daq(2)).product(5:6),'16')
          eval(['[' sprintf('AllSNs{%d},',1:(length(devices)-1)) 'AllSNs{' int2str(length(devices)) '}] = deal(devices.serialNumber);']);
          NumInterfaces = length(cell2mat(strfind(AllSNs,devices(daq(2)).serialNumber)));
          if NumInterfaces ~= 7
            ConfirmInfo(sprintf('Found %d interfaces where there should be 7.  Run "help DaqReset" for suggestions; then re-run "DaqTest"',NumInterfaces));
            return;
          end
        end
        fprintf('I''m about to try device with serial number: %s\n',devices(daq(2)).serialNumber);
        Confirm('Look again at the LED on the device you want to test and then click this button:');
        err=DaqBlinkLED(daq(2));
        if err.n
          fprintf(['\nUh oh... there was a problem when I told your second device to blink its LED.\n' ...
                   'This is the error message I received:\n\n']);

          err              
        end
        TheLEDBlinked = TwoStateQuery('Did the LED blink on the device you want to test?');
        if TheLEDBlinked < 0
          fprintf('\nYou should have clicked "Yes" or "No".  I am taking your response as a "No".\n');
          TheLEDBlinked = 0;
        end
        if TheLEDBlinked
          daq = daq(2);
          break;
        end
      end
      numTries = numTries+1;
      if numTries > 1 && ~TheLEDBlinked
        ConfirmInfo('I don''t know why this didn''t work.  I hope you do because I am bailing!');
        return;
      end
    end % while ~TheLEDBlinked
  otherwise,
    fprintf('Yay. You have %d daqs (the first one is a %s): \n',length(daq),devices(daq(1)).product);
    numTries = 0;
    TheLEDBlinked = 0;
    while ~TheLEDBlinked
      for k=1:length(daqs)
        if strcmp(devices(daq(k)).product(5:6),'16')
          eval(['[' sprintf('AllSNs{%d},',1:(length(devices)-1)) 'AllSNs{' int2str(length(devices)) '}] = deal(devices.serialNumber);']);
          NumInterfaces = length(cell2mat(strfind(AllSNs,devices(daq(k)).serialNumber)));
          if NumInterfaces ~= 7
            ConfirmInfo(sprintf('Found %d interfaces where there should be 7.  Run "help DaqReset" for suggestions; then re-run "DaqTest"',NumInterfaces));
            return;
          end
        end
        fprintf('About to test device %d, serial number: %s\n',k,devices(daq(k)).serialNumber);
        ConfirmInfo('Look at the LED on the device you want to test and click "Okay" when ready');
        err=DaqBlinkLED(daq(k));
        if err.n
          fprintf('Uh-oh... there was an error when I tried to blink the LED on device %d.  It was:\n\n',k);
          err
          fprintf('I will proceed anyways...');
        end
        TheLEDBlinked = TwoStateQuery('Did the LED blink on the device you wanted to test?');
        if TheLEDBlinked < 0
          fprintf('You should have clicked "yes" or "no".  I am taking your response as "no"');
          TheLEDBlinked = 0;
        end
        if TheLEDBlinked
          daq = daq(k);
          break;
        end
        if numTries > 1
          ConfirmInfo('I tried each LED twice...  I am giving up now.  Sorry!');
          return;
        end
        numTries = numTries+1;
      end % for k=1:length(daqs)
    end % while ~TheLEDBlinked
end % switch length(daq)

if strcmp(devices(daq).product(5:6),'16')
  Is1608=1;
else
  Is1608=0;
end

% After passing through the above code, daq should be a scalar.  References to
% it as "daq(1)" are a legacy of prior versions.  Not deemed worth fixing.  --
% mpr

if debugging
  % line numbering is... well, I think this used to be line 67...
% 	fprintf('TestClearMex once crashed at line 67.\n');
    TestClearMex(daq);
end

% fprintf(['\n** NOT RESPONDING? If your %s is not responding, try quitting and\n'...
%     '** restarting MATLAB. Denis found that this reliably restored normal\n'...
%     '** communication for his device.\n'],devices(daq(1)).product);
if strcmp(devices(daq(1)).product,'USB-1208FS') || strcmp(devices(daq(1)).product,'PMD-1208FS') || strcmp(devices(daq(1)).product,'USB-1408FS')  
  fprintf(['\n'...
      '** WIRING THE PINS OF THE USB-1208FS (and presumably -1408FS)\n'...
      '**    For the analog i/o test, please connect both channels of analog\n'...
      '** output to the first two channels of differential analog input:\n'...
      '** Connect pin 1 to 13, 2 to 12, 4 to 14, and 5 to 12.\n'...
      '**    If you''ve got one, you may want to attach a speaker or oscilloscope\n'...
      '** to analog output channel 0 (pin 13) and analog ground (pin 12). \n'...
      '**    Or you may want to connect to bit 0 of digital output port A (pin 21) \n'...
      '** and digital ground (pin 29).\n']);
elseif strcmp(devices(daq(1)).product,'USB-1608FS')
   fprintf(['\n** If you want to run all the tests we have here, you should connect analog\n' ...
            '** port 0 (pin 1) to the first bit of the digital port (pin 21).  If you would\n' ...
            '** like to test your counter, you should also connect pin 21 to the counter\n' ...
            '** terminal (pin 38).  Measurement Computing recommends that you ground any\n' ...
            '** analog ports that you are not using.  To follow that advice, connect pin 3\n' ...
            '** to pin 4, pin 5 to pin 6, pin 7 to pin 8, pin 9 to pin 10, pin 11 to pin 12,\n' ...
            '** pin 13 to pin 14, and pin 15 to pin 16.\n\n']);
else
  fprintf('Your device has not been tested with this toolbox.  You should verify the wiring and edit DaqTest to comply with it.\n');
  return;
end

ConfirmInfo(['A lot of information will be written to the Matlab command window.  ' ...
             'You can ignore most of it, but you should look at the command window now.']);

drawnow;

% The USB-1208FS is a composite device. It has four interfaces (0,1,2,3)
% which are each represented as a device in the list returned by
% PsychHID('Devices'). In calling the Daq functions, the user supplies us
% only the device index "device" corresponding to interface 0 of the
% desired USB-1208FS. However, the reports containing the samples arrive on
% interfaces 1,2,3. As returned by PsychHID('Devices'), interface i is at
% device-i, and we proceed on that basis after doing a quick check to
% confirm our assumption. However, to be platform-independent, it would be
% better to actually find all four device interfaces and confirm their
% interface numbers. USB Probe reports interface numbers, but, as far as I
% can tell, Apple's HID Explorer and our PsychHID do not return this 
% information. However, PsychHID('Devices') does report the number of outputs: 
% the USB-1208FS interface 0 has 229 (pre-Tiger) 70 (Tiger) outputs, interface 1
% has 65 (pre-Tiger) 1 (Tiger) outputs, and interfaces 2 and 3 have zero 
% outputs.  I have no idea why the number of outputs changed with the arrival of
% Mac OS X Tiger.
% 
% The USB-1608FS under Leopard behaves quite differently.  It has seven
% interfaces, no six... no five interfaces... so far I've had the system come
% up each of those ways and in that order...  In all cases, though, all but the 
% last (read highest index in the array returned by "PsychHID('Devices')") of 
% the interfaces have had 1 input, no outputs, and one collection.  The last
% interface (I am currently presuming interface 0) has 18 inputs, 66 outputs, 
% and 1 collection for a total of 85 totalElements.  I tried Tiger once and 
% found six interfaces.  Then I rebooted in Leopard again and found six 
% interfaces again as well.  Only difference seemed to be that Leopard 
% recognized the Apple Studio Display in my call to PsychHID, and Tiger didn't.
devices=PsychHID('Devices');
for d=daq 
  ok= (d-3>=1 && (devices(d-1).outputs==65 || devices(d-1).outputs==1)) || (d+3<=length(devices) && (devices(d+1).outputs==65 || devices(d+1).outputs==1));
  if devices(d).outputs<70 && ok
    fprintf(['If you have a 1408FS, you should just have been given a message\n' ...
             'suggesting some changes in the software.  If you are seeing this \n' ...
             'message, you are probably the first person to test a 1408 with the\n' ...
             'daq toolbox.\n']);
    error(sprintf('Invalid device, not the original USB-1208FS or USB-1608FS.'));
  end
end

if debugging
    TestClearMex(daq);
end

% Disable EEPROM read/write tests by default. Having this on by default is
% just bat-shit crazy!!!
if 0
  % Denis originally did not test the memory read and write functions because he
  % was afraid of messing up his firmware.  However, since this was one of the
  % first sets of functions from which I got interpretable test results, I added
  % some safeguards to DaqMemWrite and added tests of the functions here. -- mpr
  fprintf('\nTesting my ability to read and write your daq device''s EEPROM...');

  [TheOriginalData,TheErrors] = DaqMemRead(daq(1),512,59);
  MemoryReadError = 0;
  for k=1:length(TheErrors)
    if TheErrors(k).n
      MemoryReadError = 1;
      TheErrors(k)
    end
  end

  if MemoryReadError
    fprintf(['\nYou should have just been told of at least one error that occurred when I\n' ...
             'tried to read the EEPROM of your daq device.  Since it failed this test, I\n' ...
             'am aborting.  Sorry things didn''t work out...\n\n']);
    return;
  end

  MemWriteTestData = 1:59;

  fprintf('\nMemory was apparently read successfully...');
  err=DaqMemWrite(daq(1),512,MemWriteTestData);
  if err.n
    fprintf(['\nUh oh... an error occurred when I tried to write to the EEPROM of your daq\n' ...
             'device.  Since that test failed, I am aborting.  I hope that I have not screwed\n' ...
             'up the memory of your device...  In any case, I am bailing out now.  Good\n' ...
             'luck fixing things!  What I was told was:\n\n']);
    err
    return;
  end

  fprintf('\nMemory was apparently written successfully...');

  [TheWrittenData,TheErrors] = DaqMemRead(daq(1),512,59);
  MemoryReadError = 0;
  for k=1:length(TheErrors)
    if TheErrors(k).n
      MemoryReadError = 1;
      TheErrors(k)
    end
  end

  if MemoryReadError
    fprintf(['\nUh oh... I spoke too soon.  As you should just have been told, there was at\n' ...
             'least one error when I tried to read from your daq device a second time.  Odd that\n' ...
             'the device passed such a test once and then failed...  but that''s how things look\n' ...
             'to me.\n']);
    return;
  end

  if ~all(TheWrittenData(2:end) == MemWriteTestData)  
    fprintf(['\nUh oh...  I spoke too soon.  The data that I just read from the EEPROM of your\n' ...
             'daq device does not match the data that I think I just wrote to it.  I do not know\n' ...
             'why that might be.  Did you unplug it (sorry for such a lame question...).\n\n']);
    return;
  end

  fprintf('\nI just read what I wrote, and it looks like I wrote what I thought I wrote!');

  err=DaqMemWrite(daq(1),512,TheOriginalData(2:end));
  if err.n
    fprintf(['\nHmmmm...  I just tried to write your original data back to the EEPROM of your daq\n' ...
             'device.  That operation apparently failed.  I hope that I did not do any serious\n' ...
             'damage here...  It is unlikely that I did, but it is also unlikely that this function\n' ...
             'would fail here when it worked before.  Sorry things didn''t go more smoothly.  This\n' ...
             'is all I can tell you about the problem:\n\n']);
    err
    fprintf(['If you want to take matters into your own hands, the data I originally read from your\n' ...
             'daq''s EEPROM was:\n\n']);
    TheOriginalData(2:end)
    return;
  end

  fprintf(['\nEEPROM of your daq device has been restored to what it was before this test began.\n' ...
           'So far so good!\n\nMoving on to test mode changes...']);
else
  fprintf('\n\nMoving on to test mode changes...');
end

TheStatus = DaqGetStatus(daq);
if isempty(TheStatus)
  fprintf(['\nFailed to receive status report...  Look inside DaqGetStatus and run\n' ...
           'individual commands to try to figure out why this did not work.  In the\n' ...
           'mean time, I will take my leave of you.\n\n']);
  return;
else
  if isfield(TheStatus,'master')
    DaqSetSync(daq,TheStatus.master);
    NewStatus = DaqGetStatus(daq);
    if xor(TheStatus.master,NewStatus.master)
      DaqSetSync(daq,~TheStatus.master);
    else
      fprintf(['\nIt looks like I had trouble flipping the master sync status bit\n' ...
               'I don''t want to mess things up, so I return control to you.\n\n']);
      return;
    end
    FinalStatus = DaqGetStatus;
    if FinalStatus.master ~= TheStatus.master
      fprintf(['\nIt looks like I had trouble re-setting the master sync status bit\n' ...
               'after I changed it.  I am going to stop now.  I hope I haven''t caused a problem.\n\n']);
      return;
    else
      fprintf('\nSuccessfully toggled sync status bit.  Moving on...');
    end
  else
    fprintf('\nThere was a status returned, but it was missing the "master" field.  I will bail.\n\n');
    return;
  end % if isfield(TheStatus,'master'); else
  if isfield(TheStatus,'rising')
    DaqSetTrigger(daq,~TheStatus.rising);
    NewStatus = DaqGetStatus(daq);
    if xor(TheStatus.rising,NewStatus.rising)
      DaqSetTrigger(daq,TheStatus.rising);
    else
      fprintf(['\nIt looks like I had trouble flipping the trigger rising status bit\n' ...
               'I don''t want to mess things up, so I return control to you.\n\n']);
      return;
    end
    FinalStatus = DaqGetStatus;
    if FinalStatus.rising ~= TheStatus.rising
      fprintf(['\nIt looks like I had trouble re-setting the trigger rising status bit\n' ...
               'after I changed it.  I am going to stop now.  I hope I haven''t caused a problem.\n\n']);
      return;
    else
      fprintf('\nSuccessfully toggled trigger status bit.  Moving on...');
    end
  else
    fprintf('\nThere was a status returned, but it was missing the "rising" field.  I will bail.\n\n');
    return;
  end % if isfield(TheStatus,'rising'); else
  if isfield(TheStatus,'program')
    DaqPrepareDownload(daq,~TheStatus.program);
    NewStatus = DaqGetStatus(daq);
    if xor(TheStatus.program,NewStatus.program)
      DaqPrepareDownload(daq,TheStatus.program);
    else
      fprintf(['\nIt looks like I had trouble flipping the update program mode bit\n' ...
               'I don''t want to mess things up, so I return control to you.\n\n']);
      return;
    end
    FinalStatus = DaqGetStatus;
    if isempty(FinalStatus)
      fprintf(['\nIt looks like I had trouble reading the updated status from the device\n' ...
               'I don''t want to mess things up, so I return control to you.\n\n']);
      return;
    end
    if FinalStatus.program ~= TheStatus.program
      fprintf(['\nIt looks like I had trouble re-setting the update program mode bit\n' ...
               'after I changed it.  I am going to stop now.  I hope I haven''t caused a problem.\n\n']);
      return;
    else
      fprintf('\nSuccessfully toggled program status bit.  Moving on...\n');
    end
  else
    fprintf('\nThere was a status returned, but it was missing the "program" field.  I will bail.\n\n');
    return;
  end % if isfield(TheStatus,'program'); else
  
end

fprintf('\n');

fprintf('\nThe following tests assess all the %s commands, reporting any errors.\n',devices(daq(1)).product);
fprintf('\n** DIGITAL COMMANDS \n** DaqDConfigPort, DaqDIn, ');

if Is1608
  fprintf('and DaqDOut.\n\nDaqDConfigPort: configuring digital port for output.\n');
  err=DaqDConfigPort(daq,0);
  if err.n
    fprintf('Attempt to configure digital port failed.  This is all I can tell you about that:\n\n');
    err
  end
else
  % DaqDConfigPort
  fprintf('DaqDOut, and DaqGetAll.\n\nDaqDConfigPort: configuring digital ports for output.\n');
  err=DaqDConfigPort(daq(1),0,0); % configuring digital port A for output
  err=DaqDConfigPort(daq(1),1,0); % configuring digital port B for output
end

% I never had the following problem, so I don't know how helpful this is... --
% mpr

% Make sure the USB-1208FS is "attached".
if streq(err.name,'kIOReturnNotAttached')
    fprintf(['\nkIOReturnNotAttached: Mac OS error message says %s is "not attached".\n'...
        'But we know it is. We suggest that you quit and restart MATLAB. In our experience\n'...
        'that always clears up the problem.\n'],devices(daq(1)).product);
    error(sprintf('%s not attached.',devices(daq(1)).product));
end

% DaqDOut
if Is1608
  fprintf(['\nDaqDOut: random test patterns on DIO port. This is a basic functionality\n'...
      'test, showing that we can indeed read back what we wrote.\n']);
  a=Randi(255);
  err=DaqDOut(daq,a);
  fprintf('DaqDOut: %3d.\n',a);  
else % if Is1608 (DaqDOut test)
  fprintf(['\nDaqDOut: random test patterns on ports A and B. This is a basic functionality\n'...
      'test, showing that we can indeed read back what we wrote.\n']);
  a=Randi(255);
  err=DaqDOut(daq(1),0,a);
  b=Randi(255);
  err=DaqDOut(daq(1),1,b);
  fprintf('DaqDOut: %3d %3d.\n',a,b);
end % if Is1608; else (DaqDOut test)

if debugging
    TestClearMex(daq);
end

% DaqDIn
if Is1608
  fprintf('\nDaqDIn: reading digital port.\n');
  data=DaqDIn(daq);
  if isempty(data)
    disp('No data received.  That seems pathological so I will stop and let you investigate the problem.');
    return;
  else
    fprintf('DaqDIn: %3d.\n',data);
    if data ~= a
      fprintf(['What I read and what I wrote are not the same.  I will stop and let you see if\nyou can figure out why.\n' ...
               'My first guess is that you have grounded one or more of the DIO pins.  Even with the port configured\n' ...
               'for output, that seems to screw things up.  Disconnect all wires from the DIO terminals and try running\n' ...
               'DaqTest again.\n\n']);
      return;
    else
      disp('Data read matches data written.  Excellent!');
    end
  end
else % if Is1608 (DaqDIn test)
  fprintf('\nDaqDIn: reading digital ports.\n');
  data=DaqDIn(daq(1));
  if ~isempty(data)
    fprintf('DaqDIn: %3d %3d.\n',data);
    if length(data) ~= 2 || ~all(data == [a b])
      disp('Mismatch between what I read and what I wrote.  I will stop to let you investigate.');
      return;
    else
      disp('Data read matches what was written.  Excellent!');
    end
  else
    disp('No data received.  That seems pathological so I will stop and let you investigate the problem.');
    return;
  end
end % if Is1608; else (DaqDIn test)

% Bit functions
if Is1608
  fprintf('\nAbout to test ability to configure, read, and write individual bits...');
  
  fprintf('\nTwo''s complement of randomly chosen value: %d\n',bitcmp(a,8));
  TheBitValues = zeros(1,8);
  for k=1:8
    TheBitValues(k) = bitget(bitcmp(a,8),k);
  end
  DaqDWriteBit(daq,0:7,TheBitValues);
  TheResults = DaqDReadBit(daq,0:7);
  if ~all(TheResults == TheBitValues)
    disp('When I tried to write and read individual bits I came up with a different answer.  Stopping for you to investigate.');
    return;
  end
  TheWrittenValue = DaqDIn(daq);
  if TheWrittenValue == bitcmp(a,8)
    fprintf(['\nYes, the value I read back is the same as the bit complement of the previously written value!\n\n' ...
             'This test was not terribly stringent, and none of the DIO tests really check to see that your\n' ...
             'device accurately reads the inputs.  If you want to test that, you can do this:\n\n' ...
             '     DaqDConfigPort(device,1);\n\n' ...
             'and then try various combinations of wiring.  To check the value being read by all terminals at once,\n' ...
             'use:\n\n' ...
             '     DaqDReadBit(device,0:7)\n\n' ...
             '(no semi-colon at the end because you want to see the output).  Any of the channels connected to\n' ...
             'ground should read 0; others should read 1.  So, for instance, if you wire pin 21 to pin 22, pin 25\n' ...
             'to pin 26, and pin 33 to pin 34, and ran the DaqReadBit command, you should get:\n\n' ...
             '          0     1     0     1     1     1     0     1\n\n' ...
             'as an answer.  Equivalently, DaqDIn should return 186 in that case.  It is also possible to configure\n' ...
             'individual bits independently.  When I tried it I got behavior I did not fully expect (for instance,\n' ...
             'configuring a bit as output doesn''t always set the value of that bit when you write to it... as if\n' ...
             'the state of the bit when read is taken from an actual reading even though it shouldn''t be.  I did not\n' ...
             'test this too carefully, but I have no current reason to believe that there are any software errors in\n' ...
             'DaqDConfigPortBit.\n\n']);
  else
    fprintf(['Uh-oh, when I tried to read the whole port again as a single 8-bit value, I got\n' ...
             'a different answer.  Stopping for you to figure out why.\n\n']);
    return;
  end
else % if is1608 (Bit functions test)
  disp('Bitwise DIO functions have not been written for -1208 and/or -1408 daqs.  If these functions');
  disp('are important to you, you should edit them (DaqDConfigPortBit, DaqDWriteBit, DaqDReadBit) and');
  disp('then edit DaqTest so that the functions are tested on your hardware.');
  fprintf('\n\n');
end % if is1608; else (Bit functions test)

if debugging
    TestClearMex(daq);
end

% DaqGetAll
if Is1608
  fprintf('\nA working version of DaqGetAll has not been produced for the 1608, so skipping this test.\n');
else
  fprintf('\nDaqGetAll: reading all analog and digital inputs.\n');
  data=DaqGetAll(daq(1));
  if isempty(data)
    fprintf('Nothing received.\n');
  % I put this here before I set program to skip this test for 1608; if you fix
  % DaqGetAll, you will still need something like this as the 1608 has only the
  % one digital port compared to the 1208's two ports.  -- mpr
  elseif Is1608
    fprintf('data.digital: %3d.\n',data.digital);
  else
    fprintf('data.digital: %3d %3d.\n',data.digital);
  end
end
if debugging
    TestClearMex(daq);
end

% DIn/DOut speed test

fprintf(['\nDaqDOut & DaqDIn: write & read back one hundred random values as fast as\n' ...
         'possible. This test assesses reliability and latency, the time for data to \n' ...
         'go from host to device and back.\n']);
       
wrong=0;
missing=0;
if Is1608
  t=GetSecs;
  for trials=1:100
    % DaqDOut
    a=Randi(255);
    err=DaqDOut(daq,a);
    
    % Tried the following command instead of the previous command to see if 
    % the device is excessively slowed by the additional wrapping material in 
    % DaqDOut necessitated by not knowing in advance what type of Daq the user
    % has.  Found that it made *NO* difference! -- mpr
    %
    % PsychHID('SetReport',daq,2,4,uint8([0 a]));
    
    % DaqDIn
    data=DaqDIn(daq,1);
    if ~isempty(data)
      if a~=data
        wrong=wrong+1;
      end
    else
      missing=missing+1;
    end
    if err.n
      break;
    end
  end
  t=GetSecs-t;
else
  t=GetSecs;
  for trials=1:100
    % DaqDOut
    a=Randi(255);
    err=DaqDOut(daq(1),0,a);
    % DaqDIn
    data=DaqDIn(daq(1),2);
    if ~isempty(data)
      if a~=data(1)
        wrong=wrong+1;
      end
    else
      missing=missing+1;
    end
    if err.n
      break;
    end
  end
  t=GetSecs-t;
end

if missing==trials
    fprintf('Nothing received.\n');
else
    fprintf('Achieved %.0f Hz, i.e. %.0f ms cycle time. %.0f errors and %.0f missing in %.0f trials.\n',...
        trials/t,1000*t/trials,wrong,missing,trials);
end
if debugging
    TestClearMex(daq);
end

% DaqDOut
if Is1608
  fprintf('\nDaqDOut: producing fast squarewave digital output.\n');
  DaqDConfigPort(daq,0); % output
  for i=1:500
    err=DaqDOut(daq,255);
    err=DaqDOut(daq,0);
    if err.n
      break;
    end
  end
else
  fprintf('\nDaqDOut: producing fast squarewave digital output on port A.\n');
  DaqDConfigPort(daq(1),0,0); % output
  for i=1:500
    err=DaqDOut(daq(1),0,255);
    err=DaqDOut(daq(1),0,0);
    if err.n
      break;
    end
  end
end % if Is1608; else (DaqDOut test II)

if debugging
    TestClearMex(daq);
end

if Is1608
  CalibrateNow=TwoStateQuery('If you would like to calibrate analog channel 0 (recommended), wire pin 1 to pin 17 now. ',{'Okay','Skip'});
  if CalibrateNow > 0
    DaqCalibrateAIn(daq,0);
  end
  DoSquareTest = TwoStateQuery('Is pin 1 wired to pin 21?');
  if DoSquareTest > 0
    
    fprintf(['\n\nPutting square wave output through Digital Port bit 0 and reading it through\n' ...
             'Analog Input Channel 0...\n']);
           
    options.lowChannel = 0;
    options.highChannel = 0;
    options.range = 1;
    options.count = Inf;
    options.f = 100;
    
    DaqDConfigPort(daq,0);
    DaqDOut(daq,0);
    TheState=0;
    
    TheParams=DaqAInScanBegin(daq,options);

    TheInitialStart = GetSecs;
    WaitSecs(0.125);

    TheStart = GetSecs;

    for k=1:11
      AlreadyFlipped=0;
      while GetSecs-TheStart < k/4
        if ~AlreadyFlipped
          TheState = rem(TheState+1,2);
          DaqDOut(daq,TheState);
          AlreadyFlipped = 1;
        end
      end
    end
    
    WaitSecs(0.1);
    [TheData,TheParams] = DaqAInScanEnd(daq,options);

    TotalTimeElapsed = GetSecs-TheInitialStart;

    fprintf('Time elapsed: %f seconds\n',TotalTimeElapsed);
    fprintf('Number of data points: %d\n',length(TheData));
    fprintf('Expected number of data points: %d\n\n',floor(TheParams.fActual*TotalTimeElapsed));

    figure('Tag','DaqTestFig'); 
    plot((0:(size(TheData,1)-1))/TheParams.fActual-TheStart+TheInitialStart,TheData,'LineWidth',3);
    xlabel('Time (Seconds)');
    ylabel('Potential (Volts) ');
    title(sprintf('This should be a 5 V square wave with\na frequency of 2 Hz sampled at 100 Hz'));
    grid on;
    drawnow;
  end    
else % Is1608 (DaqAOut tests)
  fprintf('\n** ANALOG COMMANDS \n** DaqAIn, DaqAOut, DaqAInScan, DaqAOutScan, DaqALoadQueue\n');

  % DaqAOut: 1 kHz squarewave
  outTick=0.0005; % s
  duration=0.5; % s
  fprintf('\nDaqAOut: producing %.0f kHz squarewave output voltage for %.1f s.\n',0.5/outTick/1000,duration);
  tuntil=GetSecs;
  for i=1:duration/outTick
    tuntil=tuntil+outTick;
    WaitSecs(tuntil-GetSecs);
    err=DaqAOut(daq(1),0,mod(i,2)); % D/A 0
    err=DaqAOut(daq(1),1,mod(i,2)); % D/A 1
    if err.n
      break;
    end
  end

  % DaqAOut & DaqAIn: 4 Hz sinewave
  f=4;
  fprintf([
      '\nDaqAOut & DaqAIn: producing %.0f Hz sinewave and triangle wave output voltages,\n'...
      'on two channels, and reading them back in through the analog inputs. We plot \n'...
      'the readings in Figure 1. (It''s just noise unless you''ve installed wires\n'...
      'connecting output to input.)\n'],f);
  start=GetSecs;
  time=0:0.01:0.5;
  x=f*2*pi*time;
  sineWave=0.5+0.5*sin(x);       % Raised sine. Positive.
  triangleWave=0.5+0.5*trin(x);  % Raised triangle. Positive.
  v0=[];
  v1=[];
  for i=1:length(time)
    WaitSecs(time(i)+start-GetSecs);
    err=DaqAOut(daq(1),0,sineWave(i));
    err=DaqAOut(daq(1),1,triangleWave(i));
    if err.n
      break;
    end
    v0(i)=DaqAIn(daq(1),0,3);
    v1(i)=DaqAIn(daq(1),1,3);
  end
  if ~isempty(v0)
    figure(1);
    plot(time,v0,'r',time,v1,'g');
    ylabel('Volts');
    xlabel('Time (s)');
    title('Two channels (sine and triangle waves) put out by DaqAOut and read back in by DaqAIn.');
    text(0.15,1.9,'This graph will be meaningless noise unless you previously');
    text(0.15,1.8,'connected output to input. This test requires 4 wires ');
    text(0.15,1.7,'connecting analog output to input. Connect pin 1 to 13, ');
    text(0.15,1.6,'2 to 12, 4 to 14, and 5 to 12.');
    figure(1);
  else
    fprintf('Nothing received.\n');
  end
  if debugging
  % fprintf('TestClearMex once crashed at line 269.\n');
    TestClearMex(daq);
  end

  % DaqAInScan
  fprintf('\nDaqAInScan: sampling analog input.\n');
  options.lowChannel=[];
  options.highChannel=[];
  options.count=1000;
  options.f=2000;
  options.immediate=0;
  options.trigger=0;
  options.print=0;
  options.channel=0:0;
  options.range=3*ones(size(options.channel));
  c=length(options.channel);
  % n is sample/channel/report
  if options.immediate
    n=1;
  else
    n=31/c;
  end
  start=GetSecs;
  [data,params]=DaqAInScan(daq(1),options);
  stop=GetSecs;
  fprintf('Got %d (of %d) sample/channel on %d channels at %.0f Hz.\n',...
      size(data,1),options.count,size(data,2),params.fActual);
  if length(params.times)>1
    t=params.times;
    fprintf('Received %.1f report/s, %.1f sample/channel/s.\n',...
        (length(t)-1)/(t(end)-t(1)),(size(data,1)-n)/(t(end)-t(1)));
  end
  if debugging
  % fprintf('TestClearMex once crashed at line 302.\n');
    TestClearMex(daq);
  end

  % DaqAInScanBegin, DaqAInScanContinue, DaqAInScanEnd
  % Works well for f up to 2048/c sample/channel/s, where c is number of
  % channels. Currently (18 April 2005) unreliable at higher rates.
  c=length(options.channel);
  for log2f=11
    options.f=2^log2f/c;
    fprintf('\nDaqAInScanBegin, DaqAInScanContinue, DaqAInScanEnd: sampling analog input. \n');
    duration=1;
    options.count=duration*options.f;
    options.immediate=0;
    params=DaqAInScanBegin(daq(1),options);
    params=DaqAInScanContinue(daq(1),options);
    params=DaqAInScanContinue(daq(1),options);
    [data,params]=DaqAInScanEnd(daq(1),options);
    fprintf('Got %d (of %d) sample/channel on %d channels at %.0f Hz.\n',size(data,1),options.count,size(data,2),params.fActual);
    if length(params.times)>1
      t=params.times;
      c=size(data,2);
      % n is sample/channel/report
      if options.immediate
        n=1;
      else
        n=31/c;
      end
      fprintf('Received %.1f report/s, %.1f sample/channel/s.\n',(length(t)-1)/(t(end)-t(1)),(size(data,1)-n)/(t(end)-t(1)));
    end
  end
  if debugging
  %	fprintf('TestClearMex once crashed at line 390.\n');
    TestClearMex(daq);
  end

  % DaqAOut & DaqAInScan: 4 Hz sinewave
  % This runs well at slow sampling rates (100 Hz). Attempting to run at high
  % rates (e.g. 1000 Hz) results in losing most of the data. This isn't a
  % bug. It's asking the USB-1208FS to do two things at once (AOut and
  % AInScan) and clearly they interfere with each other. None of the
  % Measurement Computing literature that I've read mentions the issue of
  % doing two things at once, so there's no promise to live up to. This
  % simultaneous in and out is an excellent way to test the capabilities of
  % the USB-1208FS without any extra equipment, but is probably not a useful
  % model for experimental work since it limits you to low sampling rates. If
  % you need to sample in and out at the same time you may want to buy two
  % USB-1208FS units and use one for input and the other for output.
  f=1;
  fprintf([
      '\nDaqAOut & DaqAInScan: producing %.0f Hz sinewave and square wave output\n'...
      'voltages, on two channels, and reading them back in through the analog inputs. \n'...
      'We plot the readings in Figure 2. (It''s just noise unless you''ve installed\n'...
      'wires connecting output to input.)\n'],f);
  time=0:0.02:5;
  x=f*2*pi*time;
  sineWave=0.5+0.5*sin(x);       % Raised sine. Positive.
  triangleWave=0.5+0.5*trin(x);  % Raised triangle. Positive.
  squareWave=round(sineWave);
  % In
  v0=[];
  v1=[];
  options.f=100;
  options.channel=[];
  options.range=[];
  options.lowChannel=0;
  options.highChannel=0;
  options.count=round(max(time)*options.f);
  options.immediate=0;
  channel=options.lowChannel:options.highChannel;
  range=3*ones(size(channel));
  for i=channel
    DaqAIn(daq(1),i,3);% set gain range
  end
  start=GetSecs;
  params=DaqAInScanBegin(daq(1),options);
  for i=1:length(time)
    WaitSecs(time(i)+start-GetSecs);
    err=DaqAOut(daq(1),0,sineWave(i));
    err=DaqAOut(daq(1),1,squareWave(i));
    if err.n
      break;
    end
    tOut(i)=GetSecs-start;
  end
  params=DaqAInScanContinue(daq(1),options);
  [data,params]=DaqAInScanEnd(daq(1),options);
  fprintf('Got %.0f (of %.0f) sample/channel on %d channels at %.0f Hz.\n',...
    size(data,1),options.count,size(data,2),params.fActual);
  if ~isempty(data)
    t=(1:size(data,1))/params.fActual;
    v0=data(:,1);
    figure(2);
    plot(t,v0,'r');
    if size(data,2)>=2
      v1=data(:,2);
      hold on
      plot(t,v1,'g');
      hold off
      title(sprintf('%.0f Hz sine & square waves put out by DaqAOut and read back in by DaqAInScan.',f));
    else
      title(sprintf('%.0f Hz sine wave put out by DaqAOut and read back in by DaqAInScan.',f));
    end
    ylabel('Volts');
    xlabel('Time (s)');
  end
  if debugging
    TestClearMex(daq);
  end

  % DaqAOutScan. Analog output scan.
  % Works well for f up to 1024 sample/channel/s on one channel, and 256
  % sample/channel/s on two channels. Currently (28 April 2005) unreliable at
  % higher rates.
  for log2f=10
    options.f=2^log2f; % sampling frequency in Hz
    duration=0.999; % s
    fSine=10; % Hz
    options.lowChannel=0;
    options.print=0;
    outTick=1/options.f; % s
    wave=0.5+0.5*sin(fSine*2*pi*(0:outTick:duration)); % Raised sine. Positive.
    if 1
      % works for f up to 1024
      options.highChannel=0;
    else
      % works for f up to 256
      options.highChannel=1;
      wave=[wave' wave'];
    end
    options.getReports=1;
    c=1+options.highChannel-options.lowChannel;
    fprintf('\nDaqAOutScan: producing %.0f Hz sine for %.0f s sampled at %.0f Hz on %d channels.\n',...
      fSine,duration,options.f,c);
    params=DaqAOutScan(daq(1),wave,options);
    if isempty(params.start)
      params.start=nan;
    end
    if isempty(params.end)
      params.end=nan;
    end
    fprintf('Sent %.0f (out of %.0f) samples/channel on %d channels in %.0f (out of %.0f) ms.\n',...
      params.countActual,length(wave),c,1000*(params.end-params.start),1000*length(wave)/options.f);
    r=(params.numberOfReportsSent-1)/(params.end-params.start);
    fprintf('Sending speed was %.0f report/s, %.0f sample/channel/s.\n',r,r*32/c);
  end
  if debugging
  %	fprintf('TestClearMex once crashed at line 499.\n');
    TestClearMex(daq);
  end

  % DaqALoadQueue.
  fprintf('\nDaqALoadQueue: setting channel 0 to gain range 0.\n');
  err=DaqALoadQueue(daq(1),0,0);
  if debugging
    TestClearMex(daq);
  end
end % Is1608; else (DaqAOut tests)

waitsecs(0.2);

fprintf('\n** COUNTER COMMANDS \n** DaqCIn, DaqCInit\n');

% DaqCIn. Reading counter
fprintf('\nDaqCIn: reading the counter.\n');
count=DaqCIn(daq(1));
if ~isempty(count)
    fprintf('Count is %d.\n',count);
else
    fprintf('Nothing received.\n');
end

if Is1608
  IShouldTestCounter = TwoStateQuery('If you would like to test the counter, make sure pin 21 is connected to pin 38.  Then click "Yes".');
  if IShouldTestCounter > 0
    DaqCInit(daq);
    StartCount = DaqCIn(daq);
    if StartCount
      error('Counter did not reset to zero.  Sorry!');
    end
    DaqDConfigPortBit(daq,0,0);
    for k=1:150
      DaqDWriteBit(daq,0,1);
      DaqDWriteBit(daq,0,0);
    end
    EndCount = DaqCIn(daq);
    if EndCount ~= 150
      fprintf(sprintf('Counter reads %d when it should read 150!  Oops!',EndCount));
      if abs(EndCount-150) > 1
        error('In my tests, some times one count was missed, so I''ll accept that, but...');
      else
        fprintf(['\nUsually my counts were exact, but some times I missed one, particularly if I ran\n' ...
                 'the square wave test.  We''ll call this a pass, but just barely.  Moving on...\n\n']);
      end
    else
      fprintf('\nOne hundred and fifty pulses sent to the counter, 150 pulses counted.  Moving on...\n\n');
    end
  end
end % if Is1608
    

if debugging
    TestClearMex(daq);
end

% DaqCInit. Resetting counter to zero.
fprintf('\nDaqCInit: resetting counter to zero.\n');
err=DaqCInit(daq(1));

% I wrote my own versions of these above before I got down here and realized
% Denis wrote tests of these functions too...  I "if 0"d most of the code below, 
% but in the indicated places, it was like that when I got here -- mpr
%
% fprintf(['\n** MISCELLANEOUS COMMANDS \n** DaqGetStatus, DaqSetCal, DaqSetSync, DaqSetTrigger\n']);


% Denis had this test here (and had DaqBlinkLED listed in the fprintf that
% immediately precedes this), but it made more sense to me to test this up
% top and use DaqBlinkLED as a way for users to choose which device to test
% if they have more than one. I recommend using my code for this function since 
% I do more with it in my tests.  -- mpr
if 0
  % DaqBlinkLED
  fprintf('\nDaqBlinkLED: blinking the LED.\n');
  err=DaqBlinkLED(daq(1));
end

if 0 % This was commented out before I got here -- mpr
    % DaqReset. Optional reset. Make sure the USB-1208FS is "attached". In
    % theory it's nice to detect anomalous behavior and fix it, by
    % resetting the device. However, this requires a call to CLEAR PsychHID to
    % re-enumerate, which may result in a crash when you next ask PsychHID to 
    % receive reports. Don't use DaqReset. Quitting and restarting MATLAB is a 
    % reliable way of reestablishing normal communication with the USB-1208FS.
    err=PsychHID('SetReport',daq(1),2,64,uint8(64)); % Blink LED
    if streq(err.name,'kIOReturnNotAttached')
        fprintf('\nMac OS error message says %s is "not attached". But we know it is.\n',devices(daq(1)).product);
        fprintf('We will attempt to reestablish communication by re-enumerating, resetting the USB, and re-enumerating again. ');
        err=DaqReset(daq(1));
        fprintf('End of reset. Hopefully your %s is working again.\n',devices(daq(1)).product);
    end
end

if 0
  for i=[1 0]
      % DaqSetSync
      fprintf('\nDaqSetSync: configure synchronization to %d. (0 for master, 1 or 2 for slave)\n',i);
      err=DaqSetSync(daq(1),i);

      % DaqGetStatus
      fprintf('\nDaqGetStatus: reading the status bits.\n');
      status=DaqGetStatus(daq(1));
      if ~isempty(status)
          if status.master
              fprintf('master.\n');
          else
              fprintf('slave.\n');
          end
      else
          fprintf('Nothing received.\n');
      end
  end
  if debugging
      TestClearMex(daq);
  end

  if Is1608
    fprintf(['\nYou can test the calibration if you want; but I didn''t see a point\n' ...
             'to the test that was run before.  It set the calibration voltage for one\n' ...
             'second, but didn''t even tell you the pins until it was testing them.  Note\n' ...
             'the pins are different (CAL pin is pin 17 (with pin 16 being the nearest\n' ...
             'ground) in the 1608FS.  If you want to monitor the voltage across these\n' ...
             'pins, go ahead and run the test at "DaqSetCal".  I think the function\n' ...
             'should work as is, but I have not tested it.  -- mpr \n\n']);
  else
    % DaqSetCal
    fprintf(['\nDaqSetCal: setting calibration voltage output to +2.5 V for 1 s. \n'...
        'CAL is pin 16; pin 15 is ground.\n']);
    err=DaqSetCal(daq(1),1); % SetCal +2.5 V
    WaitSecs(1);
    err=DaqSetCal(daq(1),0); % SetCal 0 V

    % DaqSetTrigger
    fprintf('\nDaqSetTrigger: setting to trigger on falling edge.\n');
    err=DaqSetTrigger(daq(1),0);
  end % if Is1608; else
end % if 0
  
% Denis originally tested DaqMemRead here, but I moved the test up front because
% it does not depend upon the wiring of the terminals.  Hence it is a more
% stringent test; if it fails, it seems to me there is no reason to test
% anything else.  This is another test I recommend you stick with my code 
% instead of reverting.-- mpr
if 0
  fprintf(['\n** MEMORY COMMANDS\n'...
      '** DaqMemRead, DaqMemWrite, DaqReadCode, DaqPrepareDownload, DaqWriteCode,\n'...
      '** DaqWriteSerialNumber\n']);
  fprintf(['\nWe only test the two read commands, since the write commands could\n'...
      'potentially screw up the firmware, and we have no immediate need for\n'...
      'writing to memory.\n']);
  if debugging
      TestClearMex(daq);
  end

  % DaqMemRead
  fprintf('\nDaqMemRead: Reading memory.\n');
  data=DaqMemRead(daq(1),hex2dec('100'),62);
  if isempty(data)
      fprintf('Nothing received.\n');
  end
end

% This is another test that seems pretty pointless to me.  It might work as is,
% but if you don't put code in, how can you know what should be there?  The fact
% that *something* is received is not a terribly stringent test.  My guess is
% that if anyone wants to use this function, they'll devise their own tests of
% it, and the tests will be a good deal more meaningful than this.  I don't feel
% a need to use the function at the moment, so I'm bailing on writing that code.
%  -- mpr
if 0
  % DaqReadCode
  fprintf('\nDaqReadCode: Reading program memory.\n');
  data=DaqReadCode(daq(1),0,62);
  if isempty(data)
      fprintf('Nothing received.\n');
  end
  if debugging
      TestClearMex(daq);
  end
end

% The following was in Denis's original version; his comment is now incorrect
% because this function is tested above.  I leave this here for posterity in
% case anyone wants to revert to Denis's version. -- mpr
if 0
    % NOT TESTED. I don't need these functions, which write to memory, and
    % I don't want to risk screwing up my USB-1208FS firmware. If you decide to
    % try this, pick your addresses carefully, after consulting the table
    % of memory types: "help daqwritecode".

    % DaqMemWrite
    fprintf('\nDaqMemWrite: Writing memory.\n');
    err=DaqMemWrite(daq(1),0,1:2);

    % DaqPrepareDownload
    fprintf('\nDaqPrepareDownload: Preparing for program memory download.\n');
    err=DaqPrepareDownload(daq(1));

    % DaqWriteCode
    fprintf('\nDaqWriteCode: Writing program memory.\n');
    err=DaqWriteCode(daq(1),0,1:2);
end

% My history here was a bit odd... My recollection is that the first time I ran
% PsychHID('Devices'), the serial number reported for my USB-1608FS was 0.  I
% ran DaqWriteSerialNumber and tried to get the value set correctly.  I wasn't
% sure it worked, but later I found that the value was set (except that the
% documentation for my device has a 9-character serial number and the number is
% truncated to eight characters.  Nevertheless, I am with Denis that it is
% probably unwise to rewrite your serial number unless you absolutely know you
% need to.  So I leave this functionally commented out. -- mpr
if 0
    % WARNING: you probably do NOT want to run this command, because it
    % will over-write your factory-assigned serial number. I tested this
    % on my USB-1208FS and it works fine. Note that the new serial number
    % won't show up until you unplug and reconnect your device.

    % DaqWriteSerialNumber
    serialString='00000123';
    err=DaqWriteSerialNumber(daq(1),serialString);
    fprintf('\nDaqWriteSerial: setting the serial number to "%s".\n',serialString);
end

if Is1608
end


fprintf('\nDone.\n');

function y=trin(x)
% y=trin(x)
% Triangle wave. Has same peaks as SIN, but is linear between those peaks.
y=4*mod(x/2/pi,1);
fy=floor(y);
i=find(fy==1 | fy==2);
y(i)=2-y(i);
i=find(fy==3);
y(i)=y(i)-4;
return

function TestClearMex(daq)
% TestClearMex is being used by dgp to track down a crash of PsychHID associated
% with CLEAR MEX. 
options.printCrashers=1;
% PsychHID('ReceiveReports',daq(1),options);
% 
% % DaqDIn
% fprintf('\nDaqDIn: reading digital ports.\n');
% data=DaqDIn(daq(1));
% if ~isempty(data)
%     fprintf('DaqDIn: %3d %3d.\n',data);
% else
%     fprintf('DaqDIn: Nothing received.\n');
% end

fprintf('CLEAR MEX\n');
clear mex;

PsychHID('ReceiveReports',daq(1),options);

% DaqDIn
fprintf('DaqDIn: reading digital ports.\n');
data=DaqDIn(daq(1));
if ~isempty(data)
    fprintf('DaqDIn: %3d %3d.\n',data);
else
    fprintf('DaqDIn: Nothing received.\n');
end
return
