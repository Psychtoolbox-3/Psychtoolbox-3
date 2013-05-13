function err = LoadClut(windowPtr,clut,startEntry,bits)
% [err]=LoadClut(windowPtrOrScreenNumber,clut,[startEntry],[bits])
% 
% Load the hardware color lookup table (CLUT) of a video screen. It uses
% Screen('LoadCLUT'), as appropriate, to leave the hardware CLUT
% containing the numbers you provide in "clut", with no transformation.
% There *are* restrictions: On Microsoft Windows you can't use CLUTs for
% animation, as the operating system requires all CLUTs to contain mono-
% tonically increasing entries. Psychtoolbox currently has no way of
% detecting the resolution of your graphics cards DAC, so unless you
% explicitely provide the DAC resolution in the optional parameter 'bits',
% it will always assume a 8-Bit DAC. This assumption is safe, but it does
% not allow you to automatically take advantage of higher resolution DACs.
% Apart from that, all pixelSizes are supported. It should works with all
% graphics cards on MacOS-X, Windows and Linux. Fully supports 8-or-more-bit
% DACs. 
% 
% We *strongly* suggest that all users use Screen('LoadNormalizedGammaTable')
% in new code. Its values range between 0.0 and 1.0 and are independent
% of DAC size, the system automatically maps the range 0.0 - 1.0 to the
% range really available on your graphics card, so no need for you (or for 
% Psychtoolbox) to know the resolution of your graphics cards DAC. You will
% always automatically benefit from the highest possible resolution of your
% graphics card.
%
% 
% FUNCTION ARGUMENTS:
% 
% The err return argument is only here for backwards compatibility to
% the old Psychtoolbox. It always will be empty.
%
% "clut", the user-supplied color table, should be a clutSizex3 matrix.
% Each row in the "clut" matrix is loaded into an RGB entry in the
% hardware CLUT. The values of the matrix elements should be integers in
% the range 0 to 2^bits-1.
% 
% The maximum clut size is 256 rows, but you can pass less rows if you only
% want to change a portion of the hardware CLUT.
% 
% "startEntry" is optional and determines which hardware CLUT entry to
% load first. Entries are numbered from 0 up. The default is 0. The first
% element of "clut", i.e. clut(1), will be loaded into hardware entry
% "startEntry".
% 
% "bits" specifies how many bits you want to write to the CLUT. Typically
% it will be 8 bits, which is the default value. If you set it to
% some other value, the range of allowable entries scales accordingly.
% Thus if you use a 10-bit CLUT, then each entry should be between 0 and
% 1023, etc.
%  
% GRAPHICS CARDS WITH MORE-THAN-8-BIT DACS:
% 
% Some ATI Radeon's have 10-bit DACs. The BITS++ adapter from Cambridge
% Research Systems has 14-bit DACs.
% http://www.crsltd.com/catalog/bits++/
% 
% See also Screen subfunctions 'LoadCLUT', 'LoadNormalizedGammaTable',
% 'ReadNormalizedGammaTable'

% 4/20/06  mk     Derived it from the OS-9 PTB's LoadClut.m dated to
%                 8/24/02. Only the argument checking and parts of the
%                 online help text have been used due to the significantly
%                 different implementation in the OpenGL Psychtoolbox.

% Check the arguments
if nargin<2 || nargin>4
	error('USAGE: LoadClut(windowPtr,clut,[startEntry],[bits])');
end

if nargin<4
	bits=8;
end

if nargin<3 || isempty(startEntry)
	startEntry=0;
end

if startEntry<0 || startEntry>255
	error('startEntry %d must be in range 0 to %d',startEntry,255);
end

if max(clut(:))>2^bits-1 || min(clut(:))<0
	error(sprintf('\"clut\" values must be in range 0 to %d',2^bits-1));
end

if size(clut,1)<1 || size(clut,1)>256-startEntry
    error('Number of rows of clut must be in range 1 to 256-startEntry');
end

if size(clut,2)~=3
    error('Number of columns of clut must be 3');
end

% Setup our empty return argument:
err=[];

% Call Screens LoadCLUT and hope the best.
Screen('LoadCLUT', windowPtr, clut, startEntry, bits);

return;
