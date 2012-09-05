function kbNameResult = KbName(arg)
% kbNameResult = KbName(arg)
% 	
% 	KbName maps between KbCheck-style keyscan codes and key names.
% 	
% 	* If arg is a string designating a key label then KbName returns the 
% 	  keycode of the indicated key.  
%
% 	* If arg is a keycode, KbName returns the label of the designated key. 
%
% 	* If no argument is supplied then KbName waits one second and then 
%     calls KbCheck.  KbName then returns a cell array holding the names of
%     all keys which were down at the time of the KbCheck call. The 
%     one-second delay preceeding the call to KbCheck avoids catching the 
%     <return> keypress used to execute the KbName function. 
%
%   * If arg is 'UnifyKeyNames', KbName will switch its internal naming
%     scheme from the operating system specific scheme (which was used in
%     the old Psychtoolboxes on MacOS-9 and on Windows) to the MacOS-X
%     naming scheme, thereby allowing to use one common naming scheme for
%     all operating systems, increasing portability of scripts. It is
%     recommended to call KbName('UnifyKeyNames'); at the beginning of each
%     new experiment script.
% 	  CAUTION: This function may contain bugs. Please report them (or fix
% 	  them) if you find some.
%
%   * If arg is 'KeyNames', KbName will print out a table of all
%     keycodes->keynames mappings.
%
%   * If arg is 'KeyNamesOSX', KbName will print out a table of all
%     keycodes->keynames mappings for MacOS-X.
%
%   * If arg is 'KeyNamesOS9', KbName will print out a table of all
%     keycodes->keynames mappings for MacOS-9.
%
%   * If arg is 'KeyNamesWindows', KbName will print out a table of all
%     keycodes->keynames mappings for M$-Windows.
%
%   * If arg is 'KeyNamesLinux', KbName will print out a table of all
%     keycodes->keynames mappings for GNU/Linux, X11.
%
% 	KbName deals with keys, not characters. See KbCheck help for an 
% 	explanation of keys, characters, and keycodes.   
% 	
%   Please note that KbName always assumes a US keyboard layout. Changing
%   the keyboard layout settings in your operating system will have no
%   effect. If a keyboard with non-US layout is connected, e.g, a german
%   keyboard layout, then certain keys may not match. E.g., on a german
%   keyboard, the 'Y' key will be reported as 'Z' key and the 'Z' key will
%   be reported as 'Y' key, because these two keys are interchanged on the
%   german keyboard wrt. the US keyboard.
%
% 	There are standard character sets, but there are no standard key 
% 	names.  The convention KbName follows is to name keys with  the primary
% 	key label printed on the key.  For example, the the "]}"  key is named
% 	"]" because "]" is the primary key label and "}" is the  shifted key
% 	function.  In the case of  labels such as "5", which appears  on two
% 	keys, the name "5" designates the "5" key on the numeric keypad  and
% 	"5%" designates the QWERTY "5" key. Here, "5" specifies the primary 
% 	label of the key and the shifted label, "%" refines the specification, 
% 	distinguishing it from keypad "5".  Keys labeled with symbols not 
% 	represented in character sets are assigned names describing those
% 	symbols  or the key function, for example the space bar is named
% 	"space" and the apple  key is named "apple".  Some keyboards have
% 	identically-labelled keys distinguished 
%   only by their positions on the keyboard, for example, left and right
%   shift  keys.  Windows operating systems more recent than Windows 95 can
%   distinguish between such keys.  To name such keys, we precede the key
%   label with either  "left_" or "right_", to create the key name.  For
%   example, the left shift key  is called "left_shift".
% 	
% 	Use KbName to make your scripts more readable and portable, using key 
% 	labels instead of keycodes, which are cryptic and vary between Mac and
% 	Windows computers.  
%
% 	For example, 
% 	
% 	yesKey = KbName('return');           
% 	[a,b,keyCode] = KbCheck;
% 	if keyCode(yesKey)
% 		flushevents('keyDown');
% 		...
% 	end;
%
% OS X _ OS9 _ Windows __________________________________________________
%
%   OS X, OS 9 and Windows versions of KbCheck return different keycodes.
%   You can mostly  overcome those differences by using KbName, but with
%   some complications:
%
%   While most keynames are shared between Windows and Macintosh, not all
%   are. Some key names are used only on Windows, and other key names are
%   used only on Macintosh. For a lists of key names common to both
%   platforms and unique to each see the comments in the  body of KbName.m.
%
%   KbName will try to use a mostly shared name mapping if you add the
%   command KbName('UnifyKeyNames'); at the top of your experiment script.
%   At least the names of special keys like cursor keys, function keys and
%   such will be shared between the operating systems then.
%
%   Your computer might be able to distinguish between identically named
%   keys.  For example, left and right shift keys, or the "enter" key on
%   the keyboard and the enter key on the numeric keypad. Which of these
%   keys it can destinguish depends on the operating system. For details,
%   see comments in the body of KbName.m.
%
%   Historically, different operating systems used different keycodes
%   because they used different types of keyboards: PS/2 for Windows, ADB
%   for OS 9, and USB for OS 9, Windows, and OSX.  KbCheck on OS X returns
%   USB keycodes. 
% 	
% _________________________________________________________________________
%
% 	See also KbCheck, KbDemo, KbWait.

%   HISTORY 
%
% 	12/16/98    awi     wrote it
% 	02/12/99    dgp     cosmetic editing of comments
% 	03/19/99    dgp     added "enter" and "function" keys. Cope with hitting multiple keys.
%   02/07/02    awi     added Windows keycodes
%   02/10/02    awi     modified to return key names within a cell array in the case
%                       where no arguments are passed to KbName and it calls KbCheck.
%   02/10/02    awi     Modifed comments for changes made to Windows version. 
%   04/10/02	awi		-Cosmetic
%						-Fixed bug where "MAC2"  case was not in quotes
%						-Fixed bug where Mac loop searching for empty cells overran max index.
%   09/27/02    awi     Wrapped an index in Find() becasue Matlab will no longer index
%						matrices with logicals. 
%   06/18/03    awi     Complete rewrite featuring recursive calls to KbName, new key names, and 
%						a shared execution path for all platforms.  The key codes have been renamed
%                       to be abbreviations of standard USB keycode names.  This breaks scripts
%                       which relied on the old names but fixes the problem
%                       that KbName returned different names on Mac and Windows. See www.usb.org
%						for a table of USB HID keycodes and key names.  
%   06/24/03    awi     Added keycodes for OS 9 KbCheck  
%   06/25/03    awi     -Lowered all keycodes by 1 to match official USB table.  Previously
%						the KbCheck keycodes were the offical USB keycodes + 1 because
%						the USB table is zero-indexed whereas Matlab matrices are one-indexed.
%						As a consequence of this change KbName and KbCheck will not handle keycode
%						0, but that is ok because keycode 0 is not really used for anything.  
%						-Moved test for logical above the the test for doubles  
%						because on OS 9 logicals are also doubles.
%						-Declared the cell arrays to be persistent to conserve 
%						resources.
%	06/25/03	awi		The comment of 09/27/02 (see above) is incorrect, Matlab does still
%						support indexing with logical vectors.  Really the problem
%						was that with Matlab 6 the C mex API had changed and mxSetLogical
%						was unsupported.  KbCheck was not returning an array marked as 
%						type logical. Mex files should use mxCreateLogicalArray() or 
%						mxCreateLogicalMatrix() instead of mxSetLogical.  
% 	
%   10/12/04    awi     Cosmetic changes to comments.
%	 10/4/05    awi     Note here cosmetic changes by dgp on unknown date between 10/12/04 and 10/4/05   
%   12/31/05    mk      Transplanted the keycode table from old WinPTB KbName to OS-X KbName.
%                       This is a hack until we have PsychHID for Windows and can share the table
%                       with OS-X table.
%   21.09.06    mk      Added new command KbName('UnifyKeyNames'), which
%                       remaps many of the Windows/Linux keynames to the OS-X naming scheme for
%                       greater portability of newly written scripts.
%   27.10.06    mk      Yet another bugfix for KbName, see forum message 5247.
%   15.03.07    mk      Minimal fix to recognize uint8 type bools on Linux
%                       to fix some issue with Matlab7 on Linux.
%    6.04.09    mk      In KbName('UnifyKeyNames') Override
%                       'LockingNumLock' to become 'NumLock' on OS/X.
%   13.06.09    mk      Remove special code for Octave on all platforms except Windows.
%   21.11.09    mk      Add bugfix for failures with multiple 'Return' keys,
%                       provided by Jochen Laubrock.
%   29.09.10    mk      Add missing definition of ',<' key to Windows version.
%   20.10.11    mk      Add missing definition of 'pause' key to Windows version.

%   TO DO
%
%   -Add feature where we accept a matrix of logical row vectors.   
%
%   -Update help documentation:  Add OS X section, Explain new key names, the feature which returns
%   the table, that k=KbName;streq(KbName(KbName(k)),k) is one test of KbName,
%   that keyboards do not distinquish between left and right versions of
%   keys though the standard specifies different keycodes, the problem with the 'a' key.  

persistent kkOSX kkOS9 kkWin kkLinux kk
persistent escapedLine

%_________________________________________________________________________
% Windows part
% This is a complete hack: The code between if IsWin and the corresponding end;
% is the code of the old KbName command of the old Windows Psychtoolbox.
%
%
% These Key names are shared between Mac and Windows

if isempty(kkOSX)
	kkOSX = cell(1,256);
	kkOS9 = cell(1,256);
	kkWin = cell(1,256);
   
   % MK: Until we have a PsychHID for GNU/Linux, we query Screen for the
   % scancode to keychar mapping table:
   if IsLinux   
    kkLinux = Screen('GetMouseHelper', -3);
   else
    kkLinux = cell(1,256);
   end;
   
   % Need to manually encode sequence '\|' so Octave doesn't complain about
   % unrecognized escaped character sequence :-I
   escapedLine = char([92, 124]);
   
   % MK: Until we have a PsychHID for windows, we use this table for the
   % Windows keycode mappings. I stole it from the old WinPTB ;)
   kkWin{65} = 'a';
   kkWin{83} = 's';
   kkWin{68} = 'd';
   kkWin{70} = 'f';
   kkWin{72} = 'h';
   kkWin{71} = 'g';
   kkWin{90} = 'z';
   kkWin{88} = 'x';
   kkWin{67} = 'c';
   kkWin{86} = 'v';
   kkWin{66} = 'b';
   kkWin{81} = 'q';
   kkWin{87} = 'w';
   kkWin{69} = 'e';
   kkWin{82} = 'r';
   kkWin{89} = 'y';
   kkWin{84} = 't';
   kkWin{49} = '1!';
   kkWin{50} = '2@';
   kkWin{51} = '3#';
   kkWin{52} = '4$';  
   kkWin{53} = '5%';
   kkWin{54} = '6^';
   kkWin{187} = '=+';
   kkWin{57} = '9(';
   kkWin{55} = '7&';
   kkWin{189} = '-_';
   kkWin{56} = '8*';
   kkWin{48} = '0)';
   kkWin{221} = ']';
   kkWin{79} = 'o';
   kkWin{85} = 'u';
   kkWin{219} = '[';
   kkWin{73} = 'i';
   kkWin{80} = 'p';
   kkWin{13} = 'return';
   kkWin{76} = 'l';
   kkWin{74} = 'j';
   kkWin{222} = char(39);       % single quote
   kkWin{75} = 'k';
   kkWin{186} = ';';
   kkWin{220} = '\\';
   kkWin{188} = ',<';
   kkWin{191} = '/?';
   kkWin{78} = 'n';
   kkWin{77} = 'm';
   kkWin{190} = '.>';
   kkWin{9} = 'tab';
   kkWin{32} = 'space';
   kkWin{192} = '`';
   kkWin{46} = 'delete';
   kkWin{27} = 'esc';
   kkWin{16} = 'shift';  % Note: Windows distinguishes between left an right shift keys.
   kkWin{20} = 'capslock';
   kkWin{17} = 'control'; % Note: Windows distinguishes between left and right control keys
   kkWin{110} = '.';      
   kkWin{106} = '*';
   kkWin{107} = '+';
   kkWin{12} = 'clear';  
   kkWin{111} = '/';
   kkWin{109} = '-';
   kkWin{96} = '0';
   kkWin{97} = '1';
   kkWin{98} = '2';
   kkWin{99} = '3';
   kkWin{100} = '4';
   kkWin{101} = '5';
   kkWin{102} = '6';
   kkWin{103} = '7';
   kkWin{104} = '8';
   kkWin{105} = '9';
   kkWin{116} = 'f5';
   kkWin{117} = 'f6';
   kkWin{118} = 'f7';
   kkWin{114} = 'f3';
   kkWin{119} = 'f8';
   kkWin{120} = 'f9';
   kkWin{122} = 'f11';
   kkWin{124} = 'f13';
   kkWin{125} = 'f14';
   kkWin{121} = 'f10';
   kkWin{123} = 'f12';
   kkWin{126} = 'f15';  
   kkWin{47} = 'help';    
   kkWin{36} = 'home';
   kkWin{33} = 'pageup';
   kkWin{115} = 'f4';
   kkWin{35} = 'end';
   kkWin{113} = 'f2';
   kkWin{34} = 'pagedown';
   kkWin{112} = 'f1';
   kkWin{37} = 'left';
   kkWin{39} = 'right';
   kkWin{40} = 'down';
   kkWin{38} = 'up';
   
   % Keynames used only on Windows
   kkWin{91} = 'windows_left';
   kkWin{92} = 'windows_right';
   kkWin{93} = 'applications';
   kkWin{108} = 'seperator';
   kkWin{127} = 'f16';
   kkWin{128} = 'f17';
   kkWin{129} = 'f18';
   kkWin{130} = 'f19';
   kkWin{131} = 'f20';
   kkWin{132} = 'f21';
   kkWin{133} = 'f22';
   kkWin{134} = 'f23';
   kkWin{135} = 'f24';
   kkWin{144} = 'numlock';
   kkWin{145} = 'scrolllock';
   kkWin{246} = 'attn';
   kkWin{247} = 'crsel';
   kkWin{248} = 'exsel';
   kkWin{251} = 'play';
   kkWin{252} = 'zoom';
   kkWin{254} = 'pa1';
   kkWin{8} = 'backspace';
   kkWin{1} = 'left_mouse';
   kkWin{2} = 'right_mouse';
   kkWin{4} = 'middle_mouse';
   kkWin{45} = 'insert';
   kkWin{18} = 'alt';
   kkWin{19} = 'Pause';                            

   % Keynames used only in Windows >95
   kkWin{160} = 'left_shift';
   kkWin{161} = 'right_shift';
   kkWin{162} = 'left_control';
   kkWin{163} = 'right_control';
   kkWin{91} = 'left_menu';
   kkWin{92} = 'right_menu';
   
   
	%OS X column                                    OS 9 column             	Win Column
	
	% kk{0} = 'Undefined (no event indicated)';									%Wait until there is PsychHID 
	kkOSX{1} = 'ErrorRollOver';                              					%for Windows then use the OS X 
	kkOSX{2} = 'POSTFail';														%table.  
	kkOSX{3} = 'ErrorUndefined';
	kkOSX{4} = 'a';                                 kkOS9{1}='a';
	kkOSX{5} = 'b';                                 kkOS9{12}='b';                                                                                                       
	kkOSX{6} = 'c';                                 kkOS9{9}='c';                                         
	kkOSX{7} = 'd';                                 kkOS9{3}='d';                                          
	kkOSX{8} = 'e';                                 kkOS9{15}='e';                                          
	kkOSX{9} = 'f';                                 kkOS9{4}='f';                                          
	kkOSX{10} = 'g';                                kkOS9{6}='g';                                         
	kkOSX{11} = 'h';                                kkOS9{5}='h'; 
	kkOSX{12} = 'i';                                kkOS9{35}='i';  
	kkOSX{13} = 'j';                                kkOS9{39}='j';
	kkOSX{14} = 'k';                                kkOS9{41}='k';
	kkOSX{15} = 'l';                                kkOS9{38}='l'; 
	kkOSX{16} = 'm';                                kkOS9{47}='m'; 
	kkOSX{17} = 'n';                                kkOS9{46}='n'; 
	kkOSX{18} = 'o';                                kkOS9{32}='o';   
	kkOSX{19} = 'p';                                kkOS9{36}='p';  
	kkOSX{20} = 'q';                                kkOS9{13}='q';  
	kkOSX{21} = 'r';                                kkOS9{16}='r'; 
	kkOSX{22} = 's';                                kkOS9{2}='s';
	kkOSX{23} = 't';                                kkOS9{18}='t';
	kkOSX{24} = 'u';                                kkOS9{33}='u';
	kkOSX{25} = 'v';                                kkOS9{10}='v'; 
	kkOSX{26} = 'w';                                kkOS9{14}='w'; 
	kkOSX{27} = 'x';                                kkOS9{8}='x'; 
	kkOSX{28} = 'y';                                kkOS9{17}='y';  
	kkOSX{29} = 'z';                                kkOS9{7}='z'; 
	kkOSX{30} = '1!';                               kkOS9{19}='1!'; 
	kkOSX{31} = '2@';                               kkOS9{20}='2@';
	kkOSX{32} = '3#';                               kkOS9{21}='3#';  
	kkOSX{33} = '4$';                               kkOS9{22}='4$';
	kkOSX{34} = '5%';                               kkOS9{24}='5%'; 
	kkOSX{35} = '6^';                               kkOS9{23}='6^';  
	kkOSX{36} = '7&';                               kkOS9{27}='7&'; 
	kkOSX{37} = '8*';                               kkOS9{29}='8*';
	kkOSX{38} = '9(';                               kkOS9{26}='9(';
	kkOSX{39} = '0)';                               kkOS9{30}='0)';
	kkOSX{40} = 'Return';                           kkOS9{37}='Return';
	kkOSX{41} = 'ESCAPE';                           kkOS9{54}='ESCAPE';
	kkOSX{42} = 'DELETE';                           kkOS9{52}='DELETE';
	kkOSX{43} = 'tab';                              kkOS9{49}='tab';
	kkOSX{44} = 'space';                            kkOS9{50}='space';
	kkOSX{45} = '-_';                               kkOS9{28}='-_';
	kkOSX{46} = '=+';                               kkOS9{25}='=+'; 
	kkOSX{47} = '[{';                               kkOS9{34}='[{'; 
	kkOSX{48} = ']}';                               kkOS9{31}=']}'; 
	kkOSX{49} = escapedLine;                        kkOS9{43}=escapedLine;
	kkOSX{50} = '#-';                               
	
	kkOSX{51} = ';:';                               kkOS9{42}=';:';   
	kkOSX{52} = '''"';                              kkOS9{40}='''"';                    
	kkOSX{53} = '`~';                               kkOS9{51}='`~';
	kkOSX{54} = ',<';                               kkOS9{44}=',<'; 
	kkOSX{55} = '.>';                               kkOS9{48}='.>';
	kkOSX{56} = '/?';                               kkOS9{45}='/?';     
	kkOSX{57} = 'CapsLock';                         kkOS9{58}='CapsLock'; %FIX if other capslock
	kkOSX{58} = 'F1';                               kkOS9{123}='F1'; 
	kkOSX{59} = 'F2';                               kkOS9{121}='F2';
	kkOSX{60} = 'F3';                               kkOS9{100}='F3'; 
	kkOSX{61} = 'F4';                               kkOS9{119}='F4'; 
	kkOSX{62} = 'F5';                               kkOS9{97}='F5';
	kkOSX{63} = 'F6';                               kkOS9{98}='F6';
	kkOSX{64} = 'F7';                               kkOS9{99}='F7';   
	kkOSX{65} = 'F8';                               kkOS9{101}='F8';     
	kkOSX{66} = 'F9';                               kkOS9{102}='F9'; 
	kkOSX{67} = 'F10';                              kkOS9{110}='F10'; 
	kkOSX{68} = 'F11';                              kkOS9{104}='F11'; 
	kkOSX{69} = 'F12';                              kkOS9{112}='F12';
	kkOSX{70} = 'PrintScreen';                       
	kkOSX{71} = 'ScrollLock';                       
	kkOSX{72} = 'Pause';                            
	kkOSX{73} = 'Insert';                           
	kkOSX{74} = 'Home';                             kkOS9{116}='Home'; 
	kkOSX{75} = 'PageUp';                           kkOS9{117}='PageUp';   
	kkOSX{76} = 'DeleteForward';                    kkOS9{118}='DeleteForward';
	kkOSX{77} = 'End';                              kkOS9{120}='End';
	kkOSX{78} = 'PageDown';                         kkOS9{122}='PageDown'; 
	kkOSX{79} = 'RightArrow';                       kkOS9{125}='RightArrow'; 
	kkOSX{80} = 'LeftArrow';                        kkOS9{124}='LeftArrow';
	kkOSX{81} = 'DownArrow';                        kkOS9{126}='DownArrow';
	kkOSX{82} = 'UpArrow';                          kkOS9{127}='UpArrow';
	kkOSX{83} = 'NumLockClear';                     kkOS9{72}='NumLockClear';
	kkOSX{84} = '/';                                kkOS9{76}='/'; 
	kkOSX{85} = '*';                                kkOS9{68}='*';
	kkOSX{86} = '-';                                kkOS9{79}='-';
	kkOSX{87} = '+';                                kkOS9{70}='+';
	kkOSX{88} = 'ENTER';                            kkOS9{77}='ENTER';   
	kkOSX{89} = '1';                                kkOS9{84}='1';
	kkOSX{90} = '2';                                kkOS9{85}='2';  
	kkOSX{91} = '3';                                kkOS9{86}='3'; 
	kkOSX{92} = '4';                                kkOS9{87}='4';
	kkOSX{93} = '5';                                kkOS9{88}='5'; 
	kkOSX{94} = '6';                                kkOS9{89}='6';  
	kkOSX{95} = '7';                                kkOS9{90}='7'; 
	kkOSX{96} = '8';                                kkOS9{92}='8'; 
	kkOSX{97} = '9';                                kkOS9{93}='9';
	kkOSX{98} = '0';                                kkOS9{83}='0'; 
	kkOSX{99} = '.';                                kkOS9{66}='.';
	
	% Non-US.  
	% Typically near the Left-Shift key in AT-102 implementations.
	kkOSX{100} = ['NonUS' escapedLine];                              
	
	% Windows key for Windows 95, and ?Compose.?
	kkOSX{101} = 'Application';                     
	
	% Reserved for typical keyboard status or keyboard errors. Sent as a member of the keyboard array. Not a physical key.
	kkOSX{102} = 'Power';                                       
	kkOSX{103} = '=';                               kkOS9{82}='=';
	kkOSX{104} = 'F13';                             kkOS9{106}='F13';
	kkOSX{105} = 'F14';                             kkOS9{108}='F14';   
	kkOSX{106} = 'F15';                             kkOS9{114}='F15'; 
	kkOSX{107} = 'F16';                            
	kkOSX{108} = 'F17';                            
	kkOSX{109} = 'F18';                            
	kkOSX{110} = 'F19';                             
	kkOSX{111} = 'F20';                           
	kkOSX{112} = 'F21';                           
	kkOSX{113} = 'F22';                            
	kkOSX{114} = 'F23';                             
	kkOSX{115} = 'F24';                            
	kkOSX{116} = 'Execute';                        
	kkOSX{117} = 'Help';                            kkOS9{115}='Help';
	kkOSX{118} = 'Menu';                            
	kkOSX{119} = 'Select';                         
	kkOSX{120} = 'Stop';                            
	kkOSX{121} = 'Again';                          
	kkOSX{122} = 'Undo';                           
	kkOSX{123} = 'Cut';                          
	kkOSX{124} = 'Copy';                            
	kkOSX{125} = 'Paste';                           
	kkOSX{126} = 'Find';                            
	kkOSX{127} = 'Mute';                           
	kkOSX{128} = 'VolumeUp';                        
	kkOSX{129} = 'VolumeDown';                      
	
	%Implemented as a locking key; sent as a toggle button. Available for legacy support; however, most systems should use the non-locking version of this key.
	kkOSX{130} = 'LockingCapsLock';                  
	
	%Implemented as a locking key; sent as a toggle button. Available for legacy support; however, most systems should use the non-locking version of this key.
	kkOSX{131} = 'LockingNumLock';                  
	
	%Implemented as a locking key; sent as a toggle button. Available for legacy support; however, most systems should use the non-locking version of this key.
	kkOSX{132} = 'LockingScrollLock';               
	
	% Keypad Comma is the appropriate usage for the Brazilian keypad period (.) key. 
	%This represents the closest possible match, and system software should do the correct mapping based on the current locale setting.
	kkOSX{133} = 'Comma';                            
	
	kkOSX{134} = 'EqualSign';                       
	kkOSX{135} = 'International1';                 
	kkOSX{136} = 'International2';                  
	kkOSX{137} = 'International3';                  
	kkOSX{138} = 'International4';                  
	kkOSX{139} = 'International5';                  
	kkOSX{140} = 'International6';                  
	kkOSX{141} = 'International7';                  
	kkOSX{142} = 'International8';                  
	kkOSX{143} = 'International9';                  
	kkOSX{144} = 'LANG1';                           
	kkOSX{145} = 'LANG2';                           
	kkOSX{146} = 'LANG3';                          
	kkOSX{147} = 'LANG4';                          
	kkOSX{148} = 'LANG5';                          
	kkOSX{149} = 'LANG6';                         
	kkOSX{150} = 'LANG7';                           
	kkOSX{151} = 'LANG8';                          
	kkOSX{152} = 'LANG9';                           
	kkOSX{153} = 'AlternateErase';                 
	kkOSX{154} = 'SysReq/Attention';               
	kkOSX{155} = 'Cancel';                         
	kkOSX{156} = 'Clear';                          
	kkOSX{157} = 'Prior';                        
	kkOSX{158} = 'Return';                         
	kkOSX{159} = 'Separator';                       
	kkOSX{160} = 'Out';                          
	kkOSX{161} = 'Oper';                        
	kkOSX{162} = 'Clear/Again';                 
	kkOSX{163} = 'CrSel/Props';                  
	kkOSX{164} = 'ExSel';                           
	kkOSX{165} = 'Undefined';                       
	kkOSX{166} = 'Undefined';                       
	kkOSX{167} = 'Undefined';                       
	kkOSX{168} = 'Undefined';                       
	kkOSX{169} = 'Undefined';                       
	kkOSX{170} = 'Undefined';                       
	kkOSX{171} = 'Undefined';                       
	kkOSX{172} = 'Undefined';                      
	kkOSX{173} = 'Undefined';                       
	kkOSX{174} = 'Undefined';                      
	kkOSX{175} = 'Undefined';                       
	kkOSX{176} = '00';                              
	kkOSX{177} = '000';                            
	kkOSX{178} = 'ThousandsSeparator';             
	kkOSX{179} = 'DecimalSeparator';               
	kkOSX{180} = 'CurrencyUnit';                    
	kkOSX{181} = 'CurrencySub-unit';                
	kkOSX{182} = '(';                              
	kkOSX{183} = ')';                               
	kkOSX{184} = '{';                               
	kkOSX{185} = '}';                               
	kkOSX{186} = 'KeypadTab';                       
	kkOSX{187} = 'KeypadBackspace';                 
	kkOSX{188} = 'KeypadA';                         
	kkOSX{189} = 'KeypadB';                        
	kkOSX{190} = 'KeypadC';                         
	kkOSX{191} = 'KeypadD';                        
	kkOSX{192} = 'KeypadE';                        
	kkOSX{193} = 'KeypadF';                         
	kkOSX{194} = 'XOR';                            
	kkOSX{195} = '^';                               
	kkOSX{196} = '%';                               
	kkOSX{197} = '<';                               
	kkOSX{198} = '>';                               
	kkOSX{199} = '&';                              
	kkOSX{200} = '&&';                              
	kkOSX{201} = '|';                              
	kkOSX{202} = '||';                            
	kkOSX{203} = ':';                            
	kkOSX{204} = '#';                               
	kkOSX{205} = 'KeypadSpace';                     
	kkOSX{206} = '@';                              
	kkOSX{207} = '!';                               
	kkOSX{208} = 'MemoryStore';                     
	kkOSX{209} = 'MemoryRecall';                   
	kkOSX{210} = 'MemoryClear';                     
	kkOSX{211} = 'MemoryAdd';                      
	kkOSX{212} = 'MemorySubtract';                 
	kkOSX{213} = 'MemoryMultiply';                 
	kkOSX{214} = 'MemoryDivide';                   
	kkOSX{215} = '+/-';                             
	kkOSX{216} = 'KeypadClear';                    
	kkOSX{217} = 'KeypadClearEntry';                
	kkOSX{218} = 'KeypadBinary';                    
	kkOSX{219} = 'KeypadOctal';                     
	kkOSX{220} = 'KeypadDecimal';                  
	kkOSX{221} = 'Undefined';                       
	kkOSX{222} = 'Undefined';                       
	kkOSX{223} = 'Undefined';                       
	kkOSX{224} = 'LeftControl';                     kkOS9{60}='LeftControl';    %double entry  
	kkOSX{225} = 'LeftShift';                       kkOS9{57}='LeftShift';      %double entry
	kkOSX{226} = 'LeftAlt';                         kkOS9{59}='LeftAlt';        %double entry
	
	%Windows key for Windows 95, and ?Compose.?  Windowing environment key, examples are Microsoft Left Win key, Mac Left Apple key, Sun Left Meta key
	kkOSX{227} = 'LeftGUI';                         kkOS9{56}='LeftGUI';        %double entry
	
	kkOSX{228} = 'RightControl';                    %kkOS9{60}='RightControl'; % FIX double entry
	kkOSX{229} = 'RightShift';                      %kkOS9{57}='RightShift'; % FIX double entry
	kkOSX{230} = 'RightAlt';                        %kkOS9{59}='RightAlt';   % FIX double entry
	kkOSX{231} = 'RightGUI';                        %kkOSX{56} ='RightGUI';  % FIX double entry                                              
	kkOSX{232} = 'Undefined';                      
	kkOSX{233} = 'Undefined';                     
	kkOSX{234} = 'Undefined';                   
	kkOSX{235} = 'Undefined';                      
	kkOSX{236} = 'Undefined';                   
	kkOSX{237} = 'Undefined';                       
	kkOSX{238} = 'Undefined';                      
	kkOSX{239} = 'Undefined';                      
	kkOSX{240} = 'Undefined';                       
	kkOSX{241} = 'Undefined';                      
	kkOSX{242} = 'Undefined';                    
	kkOSX{243} = 'Undefined';                      
	kkOSX{244} = 'Undefined';                    
	kkOSX{245} = 'Undefined';                       
	kkOSX{246} = 'Undefined';                      
	kkOSX{247} = 'Undefined';                       
	kkOSX{248} = 'Undefined';                      
	kkOSX{249} = 'Undefined';                       
	kkOSX{250} = 'Undefined';                       
	kkOSX{251} = 'Undefined';                       
	kkOSX{252} = 'Undefined';                      
	kkOSX{253} = 'Undefined';                     
	kkOSX{254} = 'Undefined';                       
	kkOSX{255} = 'Undefined';                       
	kkOSX{256} = 'Undefined';                                              
	% 257-65535 E8-FFFF Reserved
	
	% Platform-specific key names.  The PowerBook G3 built-in keyboard might
	% not be 
	kkOS9{64}='MacPowerbookG3Function';
	kkOS9{53}='MacPowerbookG3Enter';
	
	% Fill in holes in the OS9 table
	for i=1:256
	    if(isempty(kkOS9{i}))
	        kkOS9{i}='Undefined';
	    end
	end

	% Choose the default table according to the platform
	if IsOS9
	        kk=kkOS9;
	elseif IsOSX
	    	kk=kkOSX;
	elseif IsWin
	        kk=kkWin;
    elseif IsLinux
            kk=kkLinux;
	end
	

end %if ~exist(kkOSX)
        
%if there are no inputs then use KbCheck to get one and call KbName on
%it.
if nargin==0
    WaitSecs(1);
    keyPressed = 0;
    while (~keyPressed)
        [keyPressed, secs, keyCodes] = KbCheck(-1); %#ok<*ASGLU>
    end
    kbNameResult= KbName(logical(keyCodes));  %note that keyCodes should be of type logical here.

elseif isempty(arg)
    % Empty argument. Could happen when the returned keyCode vector of
    % KbCheck did not report any depressed keys:
    kbNameResult=[];

%if the argument is a logical array then convert to a list of doubles and
%recur on the result. 
%Note that this case must come before the test for double below.  In Matlab 5 logicals are also
%doubles but in Matlab 6.5 logicals are not doubles.  
elseif islogical(arg) || (isa(arg,'double') && length(arg)==256) || (isa(arg,'uint8') && length(arg)==256)
    kbNameResult=KbName(find(arg));

%if the argument is a single double or a list of doubles (list of keycodes)
%or a single uint8 or list of uint8's (list of keycodes).
elseif isa(arg,'double') || isa(arg,'uint8')
    %single element, the base case, we look up the name.
    if length(arg) == 1
        if(arg < 1 || arg > 65535)
            error('Argument exceeded allowable range of 1-65536');
        elseif arg > 255 
            kbNameResult='Undefined';
        else
            kbNameResult=kk{arg};
        end;
    else
        %multiple numerical values, we iterate accross the list and recur
        %on each element.
        for i = 1:length(arg)
            kbNameResult{i}=KbName(arg(i));
        end
    end

%argument is  a single string so either it is a...
% - command requesting a table, so return the table.
% - key name, so lookup and return the corresponding key code.
elseif ischar(arg)      % argument is a character, so find the code
    if strcmpi(arg, 'Undefined')
        kbNameResult=[];            % is is not certain what we should do in this case.  It might be better to issue an error.
    elseif strcmpi(arg, 'KeyNames')  %list all keynames for this platform
        kbNameResult=kk;
    elseif strcmpi(arg, 'KeyNamesOSX')  %list all kenames for the OS X platform
        kbNameResult=kkOSX;
    elseif strcmpi(arg, 'KeyNamesOS9')
        kbNameResult=kkOS9;
    elseif strcmpi(arg, 'KeyNamesWindows')
        kbNameResult=kkWin;
    elseif strcmpi(arg, 'KeyNamesLinux')
        kbNameResult=kkLinux;
    elseif strcmpi(arg, 'UnifyKeyNames')
        % Calling code requests that we use the OS-X keyboard naming scheme
        % on all platforms. The OS-X scheme is modelled after the official
        % naming scheme for USB-HID Human interface device keyboards.
        % If we ever have a PsychHID implementation for Windows and Linux,
        % we'll have a unified keycode->keyname mapping and can get rid of
        % all this remapping cruft and the other keyboard tables.
        %
        % On OS-X and OS-9 this is a no-op, all other platforms need remapping...
        
        if IsWin
            % The following routine remaps specific Windows keycodes to their
            % corresponding OS-X / USB-HID keynames. We remap the original
            % Windows keynames as closely as possible, but there will be
            % certainly some omissions and mistakes.
            kk{8} = 'BackSpace';
            kk{13} = 'Return';
            kk{219} = '[{';
            kk{221} = ']}';
            kk{192} = '`~';
            kk{46} = 'DELETE';
            kk{27} = 'ESCAPE';
            kk{12} = 'Clear';
            kk{16} = 'Shift';
            kk{20} = 'CapsLock';
            kk{112} = 'F1';
            kk{113} = 'F2';
            kk{114} = 'F3';
            kk{115} = 'F4';
            kk{116} = 'F5';
            kk{117} = 'F6';
            kk{118} = 'F7';
            kk{119} = 'F8';
            kk{120} = 'F9';
            kk{121} = 'F10';
            kk{122} = 'F11';
            kk{123} = 'F12';
            kk{124} = 'F13';
            kk{125} = 'F14';
            kk{126} = 'F15';
            kk{160} = 'LeftShift';
            kk{161} = 'RightShift';
            kk{162} = 'LeftControl';
            kk{163} = 'RightControl';
            kk{91} = 'LeftMenu';
            kk{92} = 'RightMenu';
            kk{47} = 'Help';
            kk{36} = 'Home';
            kk{33} = 'PageUp';
            kk{45} = 'Insert';
            kk{35} = 'End';
            kk{34} = 'PageDown';
            kk{37} = 'LeftArrow';
            kk{39} = 'RightArrow';
            kk{40} = 'DownArrow';
            kk{38} = 'UpArrow';
            kk{164} = 'LeftAlt';
            kk{165} = 'RightAlt';
            kk{144} = 'NumLock';
            kk{145} = 'ScrollLock';
            kk{44} = 'PrintScreen';
            kk{91} = 'LeftGUI';
            kk{92} = 'RightGUI';
            kk{93} = 'Application';
            kk{19} = 'Pause';
        end
        
        if IsLinux
            % Remapping of Linux aka X11 keynames to OS-X/USB-HID keynames:
            % All relevant/important keys should be there now.
            KEYREMAP_TABLE = {
             'Up', 'UpArrow';
             'Down', 'DownArrow';
             'Left', 'LeftArrow';
             'Right', 'RightArrow';
             'Shift_R', 'RightShift';
             'Shift_L', 'LeftShift';
             'Prior', 'PageUp';
             'Next', 'PageDown';
             'Delete', 'DELETE';
             'Escape', 'ESCAPE';
             'Caps_Lock', 'CapsLock';
             'Control_R', 'RightControl';
             'Control_L', 'LeftControl';
             'Alt_L', 'LeftAlt';
             'Alt_R', 'RightAlt';
             'Mode_switch', 'RightAlt';
             'Super_L', 'LeftGUI';
             'Super_R', 'RightGUI';
             'Menu', 'Application';
             'Num_Lock', 'NumLock';
             'Scroll_Lock', 'ScrollLock';
             'Print', 'PrintScreen';
             'backslash', escapedLine;
             'Tab', 'tab';
             'apostrophe', '''"';
             'semicolon', ';:';
             'period', '.>';
             'comma', ',<';
             'slash', '/?';
	     'equal', '=+';
	     'minus', '-_';
	     'bracketright', ']}';
	     'bracketleft', '[{';
	     'grave', '`~';
	     'KP_Enter', 'Return';
	     'KP_Add', '+';
	     'KP_Subtract', '-';
	     'KP_Multiply', '*';
	     'KP_Divide', '/';
	     'KP_Delete', '.';
	     'KP_Insert', '0';
	     'KP_End',   '1';
	     'KP_Down',  '2';
	     'KP_Next',  '3';
	     'KP_Left',  '4';
	     'KP_Begin', '5';
	     'KP_Right', '6';
	     'KP_Home',  '7';
	     'KP_Up',    '8';
	     'KP_Prior', '9';
             '0', '0)';
             '1', '1!';
             '2', '2@';
             '3', '3#';
             '4', '4$';
             '5', '5%';
             '6', '6^';
             '7', '7&';
             '8', '8*';
             '9', '9(' };

            for i=1:length(KEYREMAP_TABLE)
                keycodes_indexes = find(strcmp(kkLinux, KEYREMAP_TABLE{i,1}));
                for index=keycodes_indexes
                    %index
                    %KEYREMAP_TABLE(i,1)
                    kk{index} = KEYREMAP_TABLE{i,2};
                end
                % If length(keycodes_indexes) == 0
                % ??? psychlasterror('reset');
                % end
            end
        end
        
        if IsOSX
            kk{131} = 'NumLock'; % Override 'LockingNumLock' ...
            % FIXME: kk{83} = 'NumLock'; % Override 'NumLockClear' as well?!?
        end
        
        % End of keyname unification code.
    else
        if IsOctave && IsWin
            % GNU/Octave on Windows does not yet support index mode for strcmpi, need to do it manually...
            for i=1:length(kk)
                if strcmpi(char(kk(i)), arg)
                    kbNameResult = i;
                    break;
                end
            end
        else
            kbNameResult=find(strcmpi(kk, arg));
        end
        if isempty(kbNameResult)
            error(['Key name "' arg '" not recognized. Maybe you need to add KbName(''UnifyKeyNames''); to top of your script?']);
        end
    end

% we have a cell arry of strings so iterate over the cell array and recur on each element.    
elseif isa(arg, 'cell')
    kbNameResult = [];
	cnt = 1;
    for i = 1:length(arg)
		codes = KbName(arg{i});
		ncodes = numel(codes);
		kbNameResult(cnt) = codes(1); %#ok<*AGROW>
		if ncodes>1
			kbNameResult(cnt+1:cnt+ncodes-1) = codes(2:ncodes);
		end
		cnt = cnt + ncodes;		
    end
else
    error('KbName can not handle the supplied argument. Please check your code or read the "help KbName".');
end
