function daq = DaqDeviceIndex(DeviceName, IShouldWarn)
% daq = DaqDeviceIndex([DeviceName][, ShowInterfaceNumberWarning=1])
% Returns a list of all your USB -1208FS, -1408FS, or -1608FS daqs.
%
% Also implements experimental code for detection of the USB-1024LS.
% However that code has not been tested yet and may need some small amount
% of tweaking to make it really work. Search the Psychtoolbox forum for
% corresponding messages...
%
% TECHNICAL NOTE: When we call PsychHID('Devices'), each USB-1208FS box
% presents itself as four HID "devices" sharing the same serial number. They
% are interfaces 0,1,2,3. They usually appear in reverse order in the
% device list. Nearly all the USB-1208FS commands use only interface 0, so
% we will select that one to represent the whole. All our Daq routines
% expect to receive just the interface 0 device as a passed designator. The
% few routines that need to access the other interfaces do so
% automatically.
%
% ADDENDUM: The above statement is correct for the 1208FS, not for the 1608FS
% and I don't know about the 1408FS.  Number of Devices found by PsychHID for a
% 1608FS and Leopard varies from five to seven with little rhyme or reason.  It
% appears that the correct number of interfaces is seven.  As with the 1208FS,
% most communication is through interface 0.  However, when acquiring data
% (e.g., using DaqAInScan), output is via interfaces 1 through 6.  Because
% PsychHID's enumeration of the interfaces is flaky, you may need to run this
% function more than once (with calls to DaqReset and probably "clear PsychHID"
% in between successive calls) in order to get device enumeration completed
% correctly.  In DaqTest, user will be warned if they try to test a 1608 and
% this function doesn't return the right number of interfaces, so I added the
% second optional argument as a flag to suppress the warning that would be
% generated here. -- mpr
%
% See also Daq, DaqFunctions, DaqPins, DaqTest, PsychHIDTest,
% DaqFind, DaqDIn, DaqDOut, DaqAIn, DaqAOut, DaqAInScan,DaqAOutScan.

% 4/15/05 dgp Wrote it.
% 8/26/05 dgp Added support for new "USB" name, as suggested by Steve Van Hooser <vanhooser@neuro.duke.edu>.
% 8/26/05 dgp Incorporated bug fixes (for compatibility with Mac OS X Tiger)
%             suggested by Jochen Laubrock <laubrock@rz.uni-potsdam.de>
%             and Maria Mckinley <parody@u.washington.edu>.
%             The reported number of outputs of the USB-1208FS has changed in Tiger.
%             http://groups.yahoo.com/group/psychtoolbox/message/3610
%             http://groups.yahoo.com/group/psychtoolbox/message/3614
%
% 11-12/xx/07  mpr Added possibility of specified input (with defaults) and
%                   changed number of expected device outputs; tested with
%                   USB-1608FS and Leopard.  Known device names: PMD-1208FS,
%                   USB-1208FS, USB-1408FS, USB-1608FS.
% 1/7/08  mpr Rewrote sections to take advantage of input, added checks to deal
%                 with problems I encountered with 1608 and PsychHID's device
%                 enumeration.
% 1/14/08 mpr added second argument
% 5/22/08 mk  Add (untested!) support for detection of USB-1024LS box.
% 5/08/13 mk  Fix enumeration of multiple connected devices (untested).
%             Use numeric USB productID instead of USB product name string
%             for matching of devices, as the latter isn't supported on
%             Windows, at least not reliably, whereas productID's are. Also
%             reject non-MCC devices based on the well-defined MCC
%             vendorID.
%             Use combined locationID and serialNumber as unique ID for
%             each instance of multiple identical DAQ devices to
%             disambiguate if user connects multiple devices of same model.
%             Bits of cleanup and refactoring.

% Flag to try to ensure that user sees warning exactly once and only if we
% have reason to believe they need to see it.
WarningFor1408 = 1;

if nargin < 2 || isempty(IShouldWarn)
    IShouldWarn=1;
end

if ~nargin || isempty(DeviceName)
    % DeviceName can be anything -- doesn't matter for further processing.
    DeviceName = '';
    AcceptAlternateNames = 1;
else
    if ~ischar(DeviceName)
        error('DaqDeviceIndex expects a character string as input.');
    end
    
    AcceptAlternateNames = 0;
    
    % Validate given DeviceName against list of supported devices:
    switch DeviceName
        case {'PMD-1208FS','USB-1208FS','PMD-1408FS','USB-1408FS'}
            if ~isempty(strfind(DeviceName,'1408'))
                fprintf(['\n\nWarning: This code has not been tested with a USB-1408FS.\n' ...
                    'You should check the number of outputs.\n\n' ...
                    'Also know that the screw terminals on the 1408FS are the \n' ...
                    'same as those on the 1208FS except that terminal 16 is \n' ...
                    'listed as "2.5 V" instead of "CAL".  Not having either of \n' ...
                    'these devices, I didn''t investigate farther.  If you test\n' ...
                    'a 1408 here, you should probably edit DaqDeviceIndex and \n' ...
                    'DaqTest indicating that the device works; fix the code if\n' ...
                    'need be to make that statement true. -- mpr\n\n']);
                % Been there; done that!
                WarningFor1408=0;
            end
        case {'PMD-1608FS','USB-1608FS'}
        case {'PMD-1024LS','USB-1024LS'}
        otherwise
            error('I did not recognize your specified device name.');
    end % switch DeviceName
end % if ~nargin || isempty(DeviceName)

NumInterfaces = [];

% Retrieve preprocessed list of PsychHID('Devices'), with some filtering
% applied to supported DAQ devices:
devices = PsychHIDDAQS;
daq=[];

for k=1:length(devices)
    % Skip all devices with a vendorId other than 2523 aka 0x09db, ie. all
    % devices that are not from Measurement Computing:
    if devices(k).vendorID ~= 2523
        continue;
    end
    
    % It is a MCC device. The product name string is only well defined on
    % OSX and Linux, but non-sensical on Windows. The USB productID is well
    % defined and accessible on all systems, so we remap the productID to
    % the name string: Database with id's available at
    % ftp://ftp.measurementcomputing.com/Customer/JBR/ULHelp/Users_Guide/Appendix/BoardType_Codes.htm
    % This remapping of product name string from productID is actually done
    % in our helper function PsychHIDDAQS above.
    
    % Set NumOutputs threshold value known to be consistent with
    % interface 0 on OSX. On Window and Linux we don't need this hack.
    %
    % Use as flag to prevent processing of unsupported devices, therefore
    % init to empty.
    NumOutputs = [];
    
    % 1024-LS type product ID?
    if devices(k).productID == 118
        % MK: This NumOutputs threshold needs to be tinkered with to find the
        % correct number for the actual output interface of 1024LS:
        NumOutputs = 0;
    end
    
    % 1608FS type product ID?
    if devices(k).productID == 125
        NumOutputs = 65;
    end
    
    % 1208FS type product ID?
    if devices(k).productID == 130
        NumOutputs = 69;
    end
    
    % 1408FS type product ID?
    if devices(k).productID == 161
        NumOutputs = 69;
    end
    
    % It is a MCC device.
    if AcceptAlternateNames
        if ~isempty(NumOutputs)
            if isempty(NumInterfaces)
                NumInterfaces = 1;
                MatchedDeviceName = devices(k).product;
                MatchedSerialNumbers = devices(k).serialNumber;
            else
                % Device already known - stored in our list?
                matchedIt = 0;
                for l=1:size(MatchedDeviceName,1)
                    if ~isempty(strfind(MatchedDeviceName(l,:),devices(k).product)) && ...
                            ~isempty(strfind(MatchedSerialNumbers(l,:),devices(k).serialNumber))
                        % Yes: Increment its interface count:
                        NumInterfaces(l) = NumInterfaces(l) + 1;
                        matchedIt = 1;
                    end
                end % for l=1:size(MatchedDeviceName,1)
                
                % Found at least 1 match in list of already known devices?
                if ~matchedIt
                    % Nope: Add this devices name/id and serialNumber to array,
                    % init to at least one interface we found by definition of
                    % reaching this point:
                    MatchedDeviceName    = char(MatchedDeviceName,devices(k).product);
                    MatchedSerialNumbers = char(MatchedSerialNumbers,devices(k).serialNumber);
                    NumInterfaces(end+1) = 1;
                end
            end % if isempty(NumInterfaces); else
            
            % We want interface #0. We get it by checking for it on Linux and
            % Windows, and indirectly detecting it by the number of outputs on OS/X:
            if (devices(k).outputs > NumOutputs) || (devices(k).interfaceID == 0)
                daq(end+1) = k;
            end
            
            if ~isempty(strfind(devices(k).product,'1408')) && WarningFor1408
                fprintf(['\n\nWarning: This code has not been tested with a USB-1408FS.\n' ...
                    'You should check the number of outputs.\n\n' ...
                    'Also know that the screw terminals on the 1408FS are the \n' ...
                    'same as those on the 1208FS except that terminal 16 is \n' ...
                    'listed as "2.5 V" instead of "CAL".  Not having either of \n' ...
                    'these devices, I didn''t investigate farther.  If you test\n' ...
                    'a 1408 here, you should probably edit DaqDeviceIndex and \n' ...
                    'DaqTest indicating that the device works; fix the code if\n' ...
                    'need be to make that statement true. -- mpr\n\n']);
                WarningFor1408=0;
            end
        end % if ~isempty(NumOutputs)
    else % if ~AcceptAlternateNames
        if streq(devices(k).product,DeviceName)
            if isempty(NumInterfaces)
                NumInterfaces = 1;
                MatchedSerialNumbers = devices(k).serialNumber;
            else
                % Device already known - stored in our list?
                matchedIt = 0;
                for l=1:size(MatchedSerialNumbers,1)
                    if ~isempty(strfind(MatchedSerialNumbers(l,:), devices(k).serialNumber))
                        % Yes: Increment its interface count:
                        NumInterfaces(l) = NumInterfaces(l) + 1;
                        matchedIt = 1;
                    end
                end
                
                % Found at least 1 match in list of already known devices?
                if ~matchedIt
                    % Nope: Add this devices serialNumber to array, init to at
                    % least one interface we found by definition of reaching
                    % this point:
                    MatchedSerialNumbers = char(MatchedSerialNumbers,devices(k).serialNumber);
                    NumInterfaces(end+1) = 1;
                end
            end % if isempty(NumInterfaces); else branch.
            
            % We want interface #0. We get it by checking for it on Linux and
            % Windows, and indirectly detecting it by the number of outputs on OS/X:
            if (devices(k).outputs > NumOutputs) || (devices(k).interfaceID == 0)
                daq(end+1) = k; %#ok<*AGROW>
            end
            
        end % if streq(devices(k).product,DeviceName)
    end % if AcceptAlternateNames; else
end % for k=1:length(devices)

if IShouldWarn
    for k=1:length(daq)
        if ~isempty(strfind(devices(daq(k)).product,'1608'))
            if NumInterfaces(k) < 7
                ConfirmInfo(sprintf(['Only found %d interfaces for DeviceIndex %d.  Execute ' ...
                    '"help DaqReset" for suggestions.'], ...
                    NumInterfaces(k),daq(k)));
            end
        end
    end
end

return;
