function PortNumber = FindSerialPort(PortString, forIOPort, dontFail)
% Syntax: PortNumber = FindSerialPort([PortString][, forIOPort=0][, dontFail=0])
%
% Purpose: Find serial port number/name associated with a particular input string.
%
% If the optional parameter 'PortString' is omitted, the routine will use a
% built-in list of default serial port names to search for. Otherwise it
% will try to find a serial port with a name matching 'PortString'.
%
% If the optional flag 'forIOPort' is zero or left out, the routine will return
% a 'PortNumber' compatible with use of the SerialComm() serial port driver
% for MacOS/X with Matlab. If 'forIOPort' is non-zero, it will return a
% serial port device string in 'PortNumber', suitable for use with the new
% IOPort serial port driver which works on all operating systems and
% runtime environments.
%
% By default, the function will abort with an error message if it doesn't
% find a matching useable port, or if it finds multiple matches, ie. if the
% mapping is ambiguous. If you set the optional flag 'dontFail' to 1, it
% will continue in such cases without any error message: It will simply
% return an empty [] PortNumber if no matching port can be found, or the
% 1st detected matching port, if multiple matches exist.
%
% Please note that support of SerialComm for OS/X is a legacy, only left
% for backwards compatibility! New code should only use the new IOPort
% driver and therefore set the 'forIOPort' flag to 1.
%
% For intersting in-depths comments on the use of this function with PR-650
% style photometers, read the source code of this file.
%

%
% History: 11/21/07   mpr   listed to port
%           1/16/08   mpr   added support for keyspan adapter
%           2/07/09   mk    Add code for use with IOPort instead of SerialComm.
%           4/05/09   mk    Add additional support for 'dontFail' = 1 flag.
%                           Add additional default 'PortString's.
%                           Add proper port detection/enumeration code for
%                           the MS-Windows operating systems. The Win32
%                           enumeration code is directly derived (with
%                           minimal modifications) from sample code donated
%                           by Xiangrui Li. Thanks!
%           6/14/09   mk    Remove redundant special case code for Octave.
%          12/04/12   zlb   Removed usa19* from PortDB since Keyserial1 is
%                           an alias for the first port anyway.
%           3/11/13   mk    Add /dev/ttyACM* devices to Linux enumeration.

% Comments from mpr:
%
% Function is slightly more general than name suggests because name reflects the
% origin of this function for finding the usb to serial port used for driving a
% device like the PR650 spectroradiometer I wrote this for.  NameOfIDString
% defaults to 'usbserial'.  I think any USB to RS-232 adapter that makes
% use of the 2303 chip (which seems to cover all of them on the market
% right now) will have the device name usbserial.  But if I'm wrong, then
% someone else will probably find out before I do.  -- mpr
%
%
% When I started writing this function I knew about PR650init, and I wanted
% to be able to use that function without knowing the port number.  I also
% found that calling that function directly leads to inconsistent
% (generally ineffective) behavior.  For some reason wrapping CMCheckInit
% solves the inconsistent/ineffective problem.  The place where I put the
% call to this function in CMCheckInit is inside a condition testing that
% user is on a machine using OS X, so the error checking here isn't
% necessary.   I also didn't know about IsOSX, IsWin, and IsOS9 when I
% wrote this -- I used the checks in DownloadPsychToolbox as my guide, and
% obviously that function can't make use of PsychToolbox functions...  So
% this could be improved even without making it work for other OSs, but I
% don't have an incentive to make improvements. -- mpr 11/21/07
%
%

if ~nargin || isempty(PortString)	
	% List of serial port devices to look for. These should all be lower case:
    if IsOSX
        PortDB = { lower('keyserial1'), lower('usbmodem'), lower('usbserial') };
    end
    
    if IsLinux
        PortDB = { lower('USB'), lower('ACM') };
    end
    
    if IsWin
        PortDB = { 'com1', 'com2', 'com3', 'com4', 'com5' };
    end
    selectionType = 'default';
else	
	% Make the port string lowercase for easier string comparison later.
	PortDB = {lower(PortString)};
	selectionType = 'requested';
end

if nargin < 2 || isempty(forIOPort)
    forIOPort = 0;
end

if nargin < 3 || isempty(dontFail)
    dontFail = 0;
end

% Init to all-empty return arguments:
PortNumber = [];

if IsWin
    % MS - Windows:
    
    % List of all candidates COM1 to COM256
    ports = strtrim(cellstr(num2str((1:256)', 'COM%i')));
    
    % For OS/X and Linux, it is easy to get all existing ports, while for
    % Windows, there seems no way to get the list. So we try all possible ports
    % one by one. Fortunately, it won't take much time if a port doesn't exist.
    
    % Start with empty cell arrays:
    availPorts={}; busyPorts={};
    
    % Disable output of IOPort during probe-run:
    oldverbosity=IOPort('Verbosity', 0);
    
    % Test each port for existence:
    for i=1:length(ports)
        
        % Try to open:
        [h, errmsg] = IOPort('OpenSerialPort', ports{i});

        % Open succeeded?
        if h>=0
            % Yes, this is an existing and available port. Close it again:
            IOPort('Close', h);
            
            % Add to list of available ports:
            availPorts{end+1} = ports{i}; %#ok<AGROW>
        elseif isempty(strfind(errmsg, 'ENOENT'))
            % Failed to open port, because it is busy, but port exists. Add
            % to list of busy ports:
            busyPorts{end+1} = ports{i}; %#ok<AGROW>
        end
        
        % Scan next COM port...
    end

    % Restore output of IOPort after probe-run:
    IOPort('Verbosity', oldverbosity);

    % Concatenate busy and ready ports into one list:
    TheTemp = [ availPorts ; busyPorts ];
    for i=1:length(TheTemp)
        ThePortDevices(i).name = char(TheTemp(i)); %#ok<AGROW>
    end

    % Done with enumeration on Windows.
else
    % Either Linux, OS/X or (hopefully) some other compatible Posix OS:

    % Get a list of any serial devices attached to the computer.
    if IsOSX
        % OS/X: Enumerate cu* devices - they correspond to serial ports:
        ThePortDevices = dir('/dev/cu*');
    else
        % Linux or some other Posix OS. Use search patterns for Linux, both
        % for USB serial converters/adaptors, and standard native serial ports:
        ThePortDevices = dir('/dev/ttyUSB*');
        ThePortDevices = [ dir('/dev/ttyACM*'); ThePortDevices];
        ThePortDevices = [ dir('/dev/ttyS*'); ThePortDevices];
    end
end

% Check to see if any serial devices were found.
if isempty(ThePortDevices) && (dontFail == 0)
	error('\n%s\n%s\n', ...
		'No serial devices were found.  Make sure the serial device is plugged', ...
		'in and the drivers installed.');
end

% For each serial type in the PortDB cell array, see if any attached serial
% devices match.
FoundIndices = [];
for i = 1:length(PortDB)
	for j = 1:length(ThePortDevices)
		if ~isempty(strfind(lower(ThePortDevices(j).name), PortDB{i}))
			FoundIndices(end+1) = j; %#ok<AGROW>
		end
	end
end
numPortsFound = length(FoundIndices);

if dontFail == 0
    if numPortsFound == 0 % No port matches.
        error('\n%s%s%s\n%s\n%s\n%s\n%s', ...
            'No ports found that match following ', selectionType, ' types(s)', ...
            makeDesiredPortTypesString(PortDB), ...
            'This is a list of the port(s) found:', ...
            makeFoundDevicesString(ThePortDevices), ...
            'Is your RS-232 adapter connected to your computer?');
    elseif numPortsFound > 1 % Too many port matches.
        error('\n%s%s%s\n%s\n\n%s\n%s', ...
            'More than one port name matched the ', selectionType, ' type(s)', ...
            makeDesiredPortTypesString(PortDB), ...
            'A more specific port type is required.', ...
            'This is a list of the port(s) found:', ...
            makeFoundDevicesString(ThePortDevices)); %#ok<SPERR>
    end
end

% Bit of a hack. If we get here and have more than 1 matching devices, we
% simply discard all but the 1st one:
if numPortsFound > 0
    % Cap to 1st entry:
    FoundIndices = FoundIndices(1);
else
    % No ports found, but not allowed to fail. We simply return, so
    % PortNumber will be an empty return argument:
    return;
end

if ~forIOPort
    % Old-Style, using SerialComm:
    % Convert the port name into a port number.

    if ~IsOSX
        error('This routine is not supported for use with SerialComm under other systems than Mac OS/X. Use with forIOPort=1 flag instead!');
    end

    if IsOctave
        error('This routine is not supported for use with SerialComm under Octave. Use with forIOPort=1 flag instead!');
    end

    if numPortsFound > 0
        portString = sprintf('/dev/%s', ThePortDevices(FoundIndices(1)).name);
        PortNumber = SerialComm('name2number', portString);
    end
else
    % New-Style: Path name of port itself for use with IOPort:
    if ~IsWin
        % Unix systems: Absolute device file path needed:
        PortNumber = sprintf('/dev/%s', ThePortDevices(FoundIndices).name);
    else
        % Windows: Pass "as is":
        PortNumber = upper(sprintf('%s', ThePortDevices(FoundIndices).name));
    end
end

return;

function s = makeDesiredPortTypesString(PortDB)
s = [];
for i = 1:length(PortDB)
	s = [char(s) sprintf('\t\t%s\n', PortDB{i})]; %#ok<AGROW>
end

function s = makeFoundDevicesString(ThePortDevices)
s = [];

if ~IsWin
    % Unices...
    if IsOSX
        prefix = 'cu.';
    else
        prefix = 'tty';
    end

    prefix = [prefix '%s'];

    for k = 1:length(ThePortDevices)
        ThisDevice = sscanf(ThePortDevices(k).name,prefix);
        s = [char(s) sprintf('\t\t%s\n',ThisDevice)]; %#ok<AGROW>
    end
else
    % MS-Windows:
    for k = 1:length(ThePortDevices)
        ThisDevice = sscanf(ThePortDevices(k).name, '%s');
        s = [char(s) sprintf('\t\t%s\n', ThisDevice)]; %#ok<AGROW>
    end    
end
