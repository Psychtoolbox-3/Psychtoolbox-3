function CMCheckInit(meterType, PortString)
% CMCheckInit([meterType], [PortString])
%
% Initialize the color meter.  The routine calls the
% lower level routines for given 'meterType'. If the low level routine
% fails, this routine retries, then prompts the user to take appropriate
% action.  If the low level routine succeeds, this
% routine  is silent.
%
% The following colormeters are supported:
%
% meterType 1 is the PR650 (default)
% meterType 2 is the CVI (need CVIToolbox) - Not yet implemented!
% meterType 3 is the CRS Colorimeter
% meterType 4 is the PR655
% meterType 5 is the PR670
%
% For the PR-series colorimeters, 'PortString' is the optional name of a
% device string for the serial port or Serial-over-USB port to which the
% device is connected. If none is given, the routine uses different files
% and built-in defaults to try to find the proper port: If there is a file
% 'CMPreferredPort.txt' in the search path of Matlab/Octave, it will parse
% that file for a PortString to use. Else it will use a hard-coded default
% inside this routine. If a calibration file with name 'PR650Ports' for the
% PR-650 colorimeter or 'PR655Ports' for the PR-655 colorimeter exists
% inside the PsychCalData folder, it will use the portname from that file.
%
% Other Colorimeters, e.g., CRS ColorCal, have their own specific setup
% methods and this routine just calls their setup code with default
% settings.
%
% 9/15/93 dhb		Modified to use new CMInit properly.
% 1/18/94 jms		Added gHardware switch
% 1/24/94 dhb		Changed sign of gHardware switch.
% 2/20/94 dhb		Call through CMETER rather than CM... routines.
% 4/4/00  dhb       Optional port name, only used on SERIAL version.
% 9/11/00 dhb       Added meterType.
% 1/4/02  dhb, ly   Try to get OS9 version to work with Megawolf board and SERIAL.
% 1/10/02 dhb, ly   Change calling convention.  Remove passing of port, but read
%                   port from a "calibration" file in PsychCalLocalData if it's there.
% 4/13/02 dgp		Cosmetic.
% 2/26/03 dhb       Add CRS Colorimeter, change meter type of PR-650 to 1.
% 10/05/06 dhb, cgb OSX version.
% 11/27/07 mpr      replaced hard coded portNameIn with FindSerialPort for OS X, and
%                   attempted to make this more robust and user-friendly.
%                   modifications tested only on Mac OS 10.5.1, Matlab 2007b, on
%                   a Mac Pro.  Other systems may require more tinkering...
% 2/07/09  mk, tbc  Add PR-655 support.
% 2/07/09  mk       Add experimental setup code for use of IOPort instead of SerialComm.
% 2/11/09  cgb,mk   Small fixes: Now we use IOPort instead of SerialComm --
%                   by default. Verified to work for both PR650 and PR655 toolboxes.
% 7/16/12  dhb      Choose right global default for 64-bit Matlab.  Hope I
%                   did this in a way that doesn't break anything else.
% 7/20/12  dhb      Undid 7/16/12 change.  Error was due to a stale IOPort on my path

global g_serialPort g_useIOPort;

% Default method of serial communication
if isempty(g_useIOPort)
    g_useIOPort = 1;
end

% Number of retries before giving up:
DefaultNumberOfTries = 5;

% Read the local preferred serial port file if it exists.
if exist(which('CMPreferredPort.txt'), 'file')
	t = textread(which('CMPreferredPort.txt'), '%s');
	defPortString = t{1};
else
    if IsOSX
        defPortString = 'usbserial';
    end
    
    if IsLinux
        defPortString = 'USB0';
    end
    
    if IsWin
        defPortString = 'COM5';
    end
end

% Set default the defaults.
switch nargin
    case 0
        meterType = 1;
        PortString = defPortString;
    case 1
        if isempty(meterType)
            meterType = 1;
        end
        PortString = defPortString;
    case 2
        if isempty(meterType)
            meterType = 1;
        end
        if isempty(PortString)
            PortString = defPortString;
        end
end

% I wrote the function FindSerialPort before I discovered CMCheckInit had
% been ported to OS X.  It may generally require less of users than relying
% on what amounts to a preference in the calibration files.  The default
% for portNameIn was 2, but on my machine it was 1 when I wrote my
% function.  Someone in the Brainard lab may want to generalize
% "FindSerialPort" to work with OSs other than Mac OS X and then use that
% function in lieu of meterports=LoadCalFile.  I am not intrepid enough to
% take that step.  -- MPR 11/21/07


switch meterType
    case 1,
        % PR-650
        % Look for port information in "calibration" file.  If
        % no special information present, then use defaults.
        meterports = LoadCalFile('PR650Ports');
        if isempty(meterports)
            if IsWin || IsOSX || IsLinux
                portNameIn = FindSerialPort(PortString, g_useIOPort);
            else
                error(['Unsupported OS ' computer]);
            end
        else
            portNameIn = meterports.in;
        end

        if IsWin || IsOSX || IsLinux
            stat = PR650init(portNameIn);
            status = sscanf(stat,'%f');
            if (isempty(status) || status == -1)
                disp('Failed initial attempt to make contact.');
                disp('If colorimeter is off, turn it on; if it is on, turn it off and then on.');
            end
            NumTries = 0;

            while (isempty(status) || status == -1) & NumTries < DefaultNumberOfTries %#ok<AND2>
                stat = PR650init(portNameIn);
                status = sscanf(stat,'%f');
                NumTries = NumTries+1;
                if (isempty(status) || status == -1) & NumTries >= 3 %#ok<AND2>
                    if IsOSX
                        if ~g_useIOPort
                            evalc(['SerialComm(''close'',' int2str(portNameIn) ');']);
                            evalc('clear SerialComm');
                        else
                            IOPort('Close', g_serialPort);
                        end
                    end
                    % Release global port handle:
                    clear global g_serialPort;
                    
                    fprintf('\n');
                    if ~rem(NumTries,4)
                        fprintf('\nHave tried making contact %d times.  Will try %d more...',NumTries,DefaultNumberOfTries-NumTries);
                    end
                end
            end
            fprintf('\n');
            if isempty(status) || status == -1
                disp('Failed to make contact.  If device is connected, try turning it off and re-trying CMCheckInit.');
            else
                disp('Successfully connected to PR-650!');
            end
        else
            error(['Unsupported OS ' computer]);
        end
    case 2,
        error('Support for CVI colormeter not yet implemented in PTB-3, sorry!');        
    case 3,
        % CRS-Colorimeter:
        if exist('CRSColorInit') %#ok<EXIST>
            CRSColorInit;
        else
            error('CRSColorInit command is missing on your path. Is the CRS color calibration toolbox set up properly?');
        end
    case 4,
        % PR-655:
        % Look for port information in "calibration" file.  If
        % no special information present, then use defaults.
        meterports = LoadCalFile('PR655Ports');
        if isempty(meterports)
            if IsWin || IsOSX || IsLinux
                portNameIn = FindSerialPort(PortString, g_useIOPort);
            else
                error(['Unsupported OS ' computer]);
            end
        else
            portNameIn = meterports.in;
        end
        
        if IsWin || IsOSX || IsLinux
            stat = PR655init(portNameIn);
            status = sscanf(stat,'%f');
            if (isempty(status) || status == -1)
                disp('Failed to make contact.  If device is connected, try turning it off and re-trying CMCheckInit.');
            else
                disp('Successfully connected to PR-655!');
            end
        else
            error(['Unsupported OS ' computer]);
		end
		
	% PR-670 - Functionality should be very similar to the PR-655, though
	%          it looks like there are a few more commands available for
	%          the 670.
	case 5
		if IsWin || IsOSX || IsLinux
			% Look for port information in "calibration" file.  If
			% no special information present, then use defaults.
			meterports = LoadCalFile('PR670Ports');
			if isempty(meterports)
				portNameIn = FindSerialPort(PortString, g_useIOPort);
			else
				portNameIn = meterports.in;
			end
			
			stat = PR670init(portNameIn);
			if strcmp(stat, ' REMOTE MODE')
				disp('Successfully connected to PR-670!');
			else
				disp('Failed to make contact.  If device is connected, try turning it off, type clear all, and then re-trying CMCheckInit.');
			end
		else
			error(['Unsupported OS ' computer]);
		end
		
    otherwise,
        error('Unknown meter type');
end
