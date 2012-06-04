% PsychSerial - Send and receive data vai serial ports.
%
% PsychSerial (formerly called SERIAL) allows you to send and receive 
% data through your computer's serial ports.
%
% Caution: This help text was written for MacOS-9, and may not contain
% accurate information about the Windows PsychSerial command, but at least
% the basic concepts should transfer...
%
% Note: PsychSerial is currently only supported on the MS-Windows
% platforms. The Linux PTB doesn't have serial port support yet and the
% MacOS-X PTB currently relies on the SerialComm driver, which has a
% slightly different syntax but basically the same functionality.
%
% Future PTB-3 releases may contain a completely rewritten PsychSerial
% driver with identical functionality on all operating system platforms,
% but for now you'll have to live with this band-aid solution.
%
% Another option is to use Matlabs 'serial' command, which however only
% works on Linux and Windows, requires the Java VM to be enabled and is
% relatively slow / high latency.
% 
% PsychSerial has multiple functions. The basic operation is to open a
% port, use it to read and write, and then close it.  For example:
% 
%   On MacOS-9:
%
%   portA = PsychSerial('Open','.Ain','.Aout',9600)
%   PsychSerial('Write',portA,['Hello world.' 13]);
%   PsychSerial('Close',portA);
%
%   On MS-Windows:
%
%   The 'outputDriverName' argument is unused, the 'inputDriverName' has a
%   different wording, e.g., 'COM1' for COM port 1 COM1.
%   portA = PsychSerial('Open','COM1','foobar',9600)
%   PsychSerial('Write',portA,['Hello world.' 13]);
%   PsychSerial('Close',portA);
%
% 
% will write the string "Hello world" followed by a carriage return to
% the A (modem) serial port at 9600 baud. 
%
% In place of a string as an argument to the 'Write' command, PsychSerial will
% accept a vector of doubles.  Each element of the vector must be an integer
% n such that 0<=n<=255.  Write interprets the doubles as ASCII values.  
% When passed a string, 'Write' will transmit all characters in the string up
% to, but not including, the first null (ascii 0).  When passed a vector of
% doubles, 'Write' will transmit all values in the vector, including nulls.
%      
% You select which port to open by specifying the names of its
% input and output drivers. 
%   ports=PsychSerial('Ports')
% will list all your ports. The most common ports are these:
%      Modem:   '.Ain', '.Aout'
%      Printer: '.Bin', '.Bout'
%
% If you add a card with additional serial ports, each port will have
% its own unique driver names. The Tips web site, below, recommends 
% suppliers of such cards, including a serial port for the Blue & 
% White G3.
% 
% web http://vision.nyu.edu/Tips/HowTo.html#UseSerial
% 
% Get a list of PsychSerial's functions by typing PsychSerial alone at the 
% Matlab prompt. Get help on any particular function by calling 
% PsychSerial with the function name followed by '?'.  For example, 
% either of these equivalent calls,
%      PsychSerial Open?
%      PsychSerial('Open?')
% provides help on the Open function.
% 
% You can use PsychSerial with multiple ports. Open ports are unavailable
% to others, so close them when you're done. PsychSerial keeps an internal 
% list of its open ports and automatically closes them all if it's 
% flushed or Matlab is terminated.
% 
% Supported baud rates: 150, 300, 600, 1200, 1800, 2400, 3600, 4800,
% 7200, 9600, 14400, 19200, 28800, 38400, 57600, 115000, 230000.
%
% You can control many of the configuration parameters, including
% two of the handshaking lines (DTR and RTS). See PsychSerial Params.
% Let us know if you need more control over handshaking.
% 
% NOTE: MacOS-9 Psychtoolbox 2.52 PsychSerial 'Read' returned all ASCII NULL 
% characters (0x00) read from the serial port.  Previously PsychSerial filtered out 
% NULLs.   This change might break some scripts which relied on the filtering.  
% Fortunately, it is easy to remove NULLs from the string returned by 
% PsychSerial, e.g., via this snippet of code:
% 
% s = PsychSerial('Read',port); 
% sWithoutNulls = s(s~=0);
%
% TROUBLESHOOTING: When making a custom cable to connect your serial
% port to another device, e.g. a photometer, it usually won't work
% until you get all the connections right. (Alas, there are certain
% lines that must be held high or low for anything to work.) Figuring
% out what's wrong is a bit tricky. The link to our Tips web site below
% will lead you to Brian Cole's quick guide to making a serial cable,
% Jim Conner's thorough explanation what all the serial port wires do,
% and let you download two serial-port debugging utilities (Serial of
% Champions and Keyboard2Serial).
% 
% web http://vision.nyu.edu/Tips/HowTo.html#UseSerial
% 
% See PsychSerialTest.m.

% 9/18/98  dhb  Wrote it.
% 8/31/99  dhb  Enhancements comment.
% 9/1/99   dgp  Added link to Tips web site.
% 9/11/99  dgp  Cosmetic.
% 10/27/99 dhb  Added note about close-order bug.
% 2/21/00  dgp  Removed close-order bug note, since it's fixed. Added high baud rates, 
% 2/21/00  dgp  Added new baud rates 115000 and 230000. 
% 4/5/00   dgp  Cosmetic. 
% 5/24/01  awi  Added paragraph explaining 'Write' with doubles vectors
% 4/28/02  awi  Changed name from "SERIAL" to "PsychSerial" to avoid conflict with Matlab's serial.  
% 1/31/08  mk   Adapted for PsychSerial on Windows with PTB-3. Not clear if
%               information here is really correct or useful though...

AssertMex('PsychSerial.m');
