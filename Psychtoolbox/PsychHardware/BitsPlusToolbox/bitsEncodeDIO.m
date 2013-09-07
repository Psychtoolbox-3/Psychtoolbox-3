% encodedDIOdata = bitsEncodeDIO(Mask, Data, Command, windowPtr, [setGammaTable])
%   Use it when you want to write DIO synchronised with screen frame.
%   It will sync timer interrupt with each frame.
%
%   'Mask' is DIO mask that must be an integer.
%
%	'Data' is a 248 element array of integers.
%
%	'Command' is the command code.
%
%	'windowPtr' is the window pointer returned by Screen('OpenWindow').
%
%	'setGammaTable' is an optional value that tells the function to set the
%	video card gamma table to a linear function if set to 1.  If set to 0,
%	the gamma table will be assumed to be linear.  By default, this is set
%	to 1.

% History:
% 18/4/05 ejw Converted it to run with OSX version of Psychtoolbox
% 12/10/2007 Modified to use common code in BitsPlusDIO2Matrix. (MK)

function encodedDIOdata = bitsEncodeDIO(Mask, Data, Command, windowPtr, setGammaTable)
% Make sure the correct number of arguments are passed.
if nargin < 4 || nargin > 5
    error('Usage: encodedDIOdata = bitsEncodeDIO(Mask, Data, Command, windowPtr, [setGammaTable])');
end

% If 'setGammaTable' isn't specified, go ahead and set it.
if nargin == 4
    setGammaTable = 1;
end

% Find out how big the window is.
[screenWidth, screenHeight] = Screen('WindowSize', windowPtr);

% Check that the screen width is at least 512 pixels wide
if screenWidth < 508
    error('window is not big enough to encode the Bits++ digital output');
end

if setGammaTable
    % THE FOLLOWING STEP IS IMPORTANT.
    % make sure the graphics card LUT is set to a linear ramp
    % (else the encoded data will not be recognised by Bits++).
    % There is a bug in some of the OpenGL drivers such that the
    % gamma is incorrectly mapped from [0, 255/256] instead of [0, 1].
    % If you are having trouble using the DIO, try uncommenting the 2nd line
    % below.
    LoadIdentityClut(windowPtr);
end

% Call common encoder routine to create T-Lock image matrix:
encodedDIOdata = BitsPlusDIO2Matrix(Mask, Data, Command);

% Display the Tlock packet on the back buffer starting at the bottom
% left.
encodedLineRect = [0, screenHeight - 1, size(encodedDIOdata, 2), screenHeight];
Screen('PutImage', windowPtr, encodedDIOdata, encodedLineRect);

% Make the make buffer current during the next blanking period
Screen('Flip', windowPtr);
