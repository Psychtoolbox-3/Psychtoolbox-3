function daqdevs = PsychHIDDAQS
% daqdevs = PsychHIDDAQS -- Enumerate and preprocess HID-DAQ device list.
%
% Used as helper by Daq toolbox functions. Retrieves the USB-HID devices
% list from PsychHID, prefilters all DAQ device entries, so they have a
% well defined productID, product name string and a somewhat unique
% serialNumber -- as far as this is possible. Leaves unsupported (non-DAQ)
% device entries alone.
%
% It then caches the preprocessed list to save processing time on
% successive calls to the function.
%

% History:
% 8-May-2013  mk  Wrote it.

% Cache devices in 'devices':
persistent devices;

% Already detected and cached?
if isempty(devices)
    % No. Load PsychHID:
    LoadPsychHID;
    
    % Retrieve device list from it:
    devices = PsychHID('Devices');

    % Iterate through device list and preprocess the MCC-DAQ devices, while
    % leaving other HID devices unmodified:
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
        
        % 1024-LS type product ID?
        if devices(k).productID == 118
            devices(k).product = 'USB-1024LS';
        end
        
        % 1608FS type product ID?
        if devices(k).productID == 125
            devices(k).product = 'USB-1608FS';
        end
        
        % 1208FS type product ID?
        if devices(k).productID == 130
            devices(k).product = 'USB-1208FS';
        end
        
        % 1408FS type product ID?
        if devices(k).productID == 161
            devices(k).product = 'USB-1408FS';
        end
        
        % Abuse serialNumber as unique Id of device: We combine the device
        % serialNumber string and its numeric USB locationID into a single
        % string to distinguish between multiple different devices of the same
        % productID. Why the combination? locationID is a relatively unique
        % value for a specific location on the USB bus for OSX and Linux,
        % probably unique unless somebody manages to connect two identical DAQ
        % devices to the same USB host port via some weird intermediate USB
        % hub. Even then it might be unique, but this is unconfirmed. However,
        % the locationID is mostly useless on Windows, despite claims otherwise
        % in MS documentation, because it has the same value for identical
        % products, even if they are connected to different host USB ports!
        % serialNumber is in theory unique for each device, but in practice
        % many vendors are lazy and don't assign a unique serialNumber, or even
        % a serialNumber at all! Net result is that both values serve as unique
        % device id sometimes, but not always. The best we can do is combine
        % both and hope that not both fail at the same time.
        devices(k).serialNumber = [devices(k).serialNumber '::' sprintf('%f', devices(k).locationID)];
    end
end

% Return cached, preprocessed device list:
daqdevs = devices;

return;
