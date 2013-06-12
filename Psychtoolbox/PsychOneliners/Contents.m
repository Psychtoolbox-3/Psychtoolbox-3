% Psychtoolbox:PsychOneliners.
% 
% help Psychtoolbox % For an overview, triple-click me & hit enter.
% 
%   AddStructs              - Merges input structs into one struct.
%   AltSize                 - ALTSIZE is an extension of SIZE, it supports querying the size of multiple dimensions of a variable in one call.
%   AreStructsEqualOnFields - Are two structures the same on the passed list of fields?
%   Ask                     - Display message, get user's response.
%   AssertGLSL              - Require the OpenGL shading language is supported.
%   AssertMex               - Detect missing mex file and exit with error.
%   AssertOpenGL            - Require that Psychtoolbox be based on OpenGL.
%   AssertOSX               - Require that Psychtoolbox be based on OS X.
%   BackupCluts             - Make internal backup of (given) Cluts, for later restoration via RestoreCluts.
%   BlackIndex              - Returns number that will produce the color black.
%   CatStr                  - Concatenate array or cell array of strings.
%   CenterMatOnPoint        - Returns indices to center a matrix on a point.
%   Circle                  - Returns a boolean mask in the shape of a Circle.
%   CleanStruct             - Deletes all empty structs and fields from a struct array, optionally recursively.
%   CreateUniformDotsIn3DFrustum - Sample dots in 3D frustum uniformly.
%   CropBlackEdges          - Detects if there are any black edges around an image and returns indices that can be used to cut away these edges.
%   DeEmptify               - Deletes empty cells or rows from cellarray.
%   DegToMrad               - Convert angle in degrees to milliradians (mrad).
%   DescribeComputer        - Print a description of the environment.
%   DotOffset               - Calculate offsets for a 3D movement. Various per-dot options.
%   Ellipse                 - Returns a boolean mask in the shape of a (super-) Ellipse.
%   EnforcePos              - Truncate negative values of a vector to 0.
%   Explode                 - Splits a numeric or character array by a delimiter or delimiter-pattern.
%   FillEmptyFields         - Fill all empty fields of a struct(array) or all empty elements of a cellarray with specified value.
%   FindInd                 - Returns indices in all dimensions to non-zero elements in matrix.
%   FindRepeatAlongDims     - Find repeated rows or columns in a matrix.
%   FunctionFolder          - Get full path to folder containing passed function.
%   GetEchoNumber           - Get a number typed on-screen.
%   GetEchoString           - Get a string typed on-screen.
%   GetKbChar               - Simple, limited replacement for GetChar(), using KbCheck for character input.
%   GetMyCaller             - Returns the name of the calling function.
%   GetNumber               - Get a number typed at the keyboard.
%   GetString               - Get a string typed at the keyboard.
%   GetSubversionPath       - Return path required to invoke snv.
%   GetSVNInfo              - Find info on SVN version number of directory or file.
%   GetWithDefault          - Get string or number with prompt and default value.
%   GrayIndex               - Any graylevel from black (0) to white (1).
%   GroupStructArrayByFields - An aid to sorting data kept in structure arrays.
%   hexstr                  - Hex string of lowest 32 bits of any number.
%   ImageToVec              - Convert a grayscale image to vector format.
%   Ind2Str                 - Converts numbers to characters (decimal to base 26 conversion). Useful for character indices.
%   Interleave              - Interleaves any number of arrays. Can handle different data types.
%   IsACell                 - Tests (recursively--cells in cells) if a cell satisfies a user-supplied condition.
%   IsARM                   - Return if running on a processor with ARM architecture, typically a mobile or embedded system.
%   IsGLES                  - Return if the current rendering api in use is OpenGL-ES, the "OpenGL Embedded Subset".
%   IsGLES1                 - Return if the current rendering api in use is OpenGL-ES 1.x.
%   IsGUI                   - Is the Matlab or Octave GUI enabled in this session?
%   IsLinux                 - Shorthand for testing whether running under Linux.
%   IsMinimumOSXVersion     - Query if this is a specific OS/X version or higher.
%   IsOctave                - Shortand for testing whether running under Octave.
%   IsOSX                   - Return if running on a Apple OSX operating system.
%   IsWin                   - Return if running on a MS-Windows operating system.
%   Is64Bit                 - Return if script is running on a 64-Bit Octave or Matlab.
%   KbMapKey                - Checks if any of specified keys is depressed in a vector returned by KbCheck, KbWait etc.
%   kPsychGUIWindow         - Flag to ask Screen() to create onscreen windows with behaviour similar to normal GUI windows.
%   LoadIdentityClut        - Loads the identity CLUT on a specified monitor.
%   log10nw                 - Compute log base 10 without annoying warnings.
%   MacModelName            - Mac model name, e.g. 'PowerBook G4 15"'.
%   Magnify2DMatrix         - Expand the size of a two-dimensional matrix via entry replication.
%   MakeBeep                - Compute a beep of specified frequency and duration, for Snd.
%   MakeCosImage            - Make a cosinusoidal image.
%   MakeGrid                - Makes raster of elements centered on screen / in image (leftover space is divided equally over the edges).
%   MakeSincImage           - Make a sinc image.
%   MakeSineImage           - Make a sinusoidal image.
%   MapIndexColorThroughClut - Convert an index color image and clut to a true color image.
%   MergeCell               - Concatenates contents of input cells element-wise.
%   MradToDeg               - Convert angle in milliradians (mrad) to degrees.
%   NameBytes               - Nicely format memory size for human readers.
%   NameFrequency           - Nicely format clock rate for human readers.
%   NearestResolution       - Find a screen resolution that most closely matches a requested resolution.
%   OSName                  - Convential English-language name of your operating system.
%   PackColorImage          - Pack three color planes into one m by n by three matrix.
%   ProgressBar             - Displays a progress bar in MATLAB's command window.
%   PsychDebugWindowConfiguration - Enable special debug window configuration to aid single display debugging.
%   PsychGPUControl         - Control low-level operating parameters of certain supported GPU's.
%   PsychNumel              - Drop-in replacement for numel() on old Matlab versions that don't support it.
%   PsychtoolboxRoot        - Robust way to get path to Psychtoolbox folder, even if renamed.
%   RemoveMatchingPaths     - Removes folders that contain a given string from the path.
%   RemoveSVNPaths          - Removes ".svn" folders from the path.
%   Replace                 - Perform exact Replace on strings or numeric arrays.
%   Resolute                - Cuts from and adds to a matrix to make it the specified dimensions.
%   RestoreCluts            - Restore original CLUT's for all monitors from backups made during LoadIdentityClut().
%   Rot3d                   - Rotates a matrix in 3D space (around the x, y or z axis) in 90 degrees steps.
%   SaveIdentityClut        - Store current or given CLUT as identity LUT for use with LoadIdentityClut.
%   SaveMovieFrames         - Displays a GUI in which a movie can be played and from which screenshots can be saved.
%   sca                     - Shorthand for Screen('CloseAll').  Using this is a good way to make your code obscure.
%   ScreenDacBits           - What is precision of the graphics boardDACs. Currently returns 8 always.
%   SetResolution           - Change display resolution, refresh rate and color depths to requested values.
%   ShrinkMatrix            - Shrinks a 2-D or 3-D matrix (an image) by a factor.
%   SmartVec                - Creates a vector/sequence that satisfies certain conditions.
%   SortCell                - Sorts cell matrices containing different data types in different columns.
%   Speek                   - Use speech synthesis output to speak a given text. Mac OS/X only.
%   Stopwatch               - Time intervals.
%   streq                   - strcmp.
%   StrPad                  - Makes a string a specified length, either by pre-padding it with a specified character or cutting from its beginning.
%   Struct2Vect             - Returns (cell-) array with all values in a specified field of a structure array.
%   TextBounds              - Draw string, return enclosing rect.
%   TextCenteredBounds      - Draw string, centered, return enclosing rect.
%   UnpackColorImage        - Extract three color planes from an m by n by 3 color image.
%   VecToImage              - Convert a grayscale image from vector to image format.
%   WhiteIndex              - Returns number that will produce the color white.
%   WinDesk                 - Sends command to windows shell to minimize all Windows, equal to Windows+M.
%   WrapString              - Word wrap (break into lines).
%
%
% The following is a list of old one-liners that might some day be updated
% from PTB-2, but haven't been yet.
%
%   BlankingInterruptRate - Used by PsychBasic FrameRate.
%   ClutDefault           - Returns standard clut for screen at any pixelSize.
%   CmdWinToUpperLeft     - Bring Command window forward, saving Screen window.
%   DescribeScreen        - Print a description of the screen's video driver.
%   DescribeScreenPrefs   - Print more about the screen's video driver.
%   GammaIdentity         - Returns an identity gamma table appropriate to the screen's dacSize.
%   IsDownArrow           - Is char the down arrow?
%   IsLeftArrow           - Is char the left arrow?
%   IsRightArrow          - Is char the right arrow?
%   IsUpArrow             - Is char the up arrow?
%   IsInOrder             - Are the two strings in alphabetical order?
%   IsPopCharProInstalled - Is the Control Panel "PopChar Pro" installed?
%   MaxPriorityGetSecs    - Figure out the maximum priority compatible with GetSecs. Slow.
%   ScreenClutSize        - How many entries in the graphic card Color Lookup Table?
%   ScreenUsesHighGammaBits - Does this card use the high 10 bits of the gamma values?
%   SCREENWinToFront      - Bring Screen window back in front of Command window.
%   ShowTiff              - Show a TIFF file, calibrated.
