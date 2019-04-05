function versionString=AppleVersion(gestaltString)
% versionString=AppleVersion(gestaltString)
%
% WARNING: This function is deprecated and will likely cease to work
% in a future Psychtoolbox release. This is due to Apple deprecating
% their Gestalt() function from their operating system. While Gestalt,
% and thereby this function, still works on OSX 10.10, no guarantees
% can be made about future OSX versions.
%
% Adapt your code accordingly to do without this function.
%
% OS 9 and OS X: __________________________________________________________
%
% AppleVersion('qtim') % QuickTime
% AppleVersion('atkv') % AppleTalk
% Uses Gestalt to retrieve Apple version information and returns a string.
% Apple has a "standard" format for versions, e.g. 3.0f7 or 8.0.1, that
% they use for some of their software components. APPLEVERSION uses GESTALT
% to retrieve the information. However, this is only useful for the few
% Gestalt selectors that return information in this "standard" format.
% If the selector is undefined (possibly because that software is
% not present) APPLEVERSION returns an empty string.
%
% WINDOWS: ________________________________________________________________
% 
% AppleVersion does not exist in Windows.
%
% _________________________________________________________________________
%
% see also: Gestalt, Screen('Computer?'), MacModelName

% Matlab bugs fixed in Matlab 5.2.1:
% eval('b=gestalt(gestaltString);','return;');
% error message is not suppressed. RETURN is ignored.
% eval('gestalt(''atkv'');','[]')
% semicolon causes first arg to be treated as error.

% AppleVersion('apvr') % gestaltAppearanceVersion
% AppleVersion('ascv') % gestaltAppleScriptVersion
% AppleVersion('csvr') % gestaltControlStripVersion
% AppleVersion('otrv') % gestaltOpenTptRemoteAccessVersion
% AppleVersion('sysu') % gestaltSystemUpdateVersion
% AppleVersion('cltn') % gestaltCollectionMgrVersion
% AppleVersion('walk') % gestaltALMVers
% AppleVersion('gestaltALMVers') % gestaltAutoBuildVersion
% #define gestaltATalkVersion 'atkv' /* AppleTalk version &AD01/M01 */
% WARNING:
% This selector returns the majorRev field of the NumVersion record as
% hexadecimal instead of the usual BCD.
% 
% #define gestaltGestaltKaputVersion 'G\0xa0K\0xa0' /* Gestalt Kaput */
% NOTE: Both the t characters are actually the option-t character
% (0xA0).
% 
% #define gestaltGestaltVersion 'G\0x8ast' /* Gestalt version */
% NOTE: The "a" is actually the option-u/a character (0x8A).
% #define gestaltUniversalDiskFormatVersion? 'kudf'


% HISTORY:
% 3/21/98   dgp wrote it.  
% 5/19/99   dgp Check for "unavailable selector" error.
% 12/7/04   awi Capitalized "Gestalt". Divided help by platform.
% 1/29/05   dgp Cosmetic.


persistent firstTime

if isempty(firstTime)
  firstTime = 0;
  warning('The Psychtoolbox function AppleVersion() is deprecated due to Apple''s fault. The function may cease to work in the future! Adapt your code accordingly.');
end

b=eval('Gestalt(gestaltString)','[]');
if isempty(b) || (isscalar(b) && b==-5551)
	versionString='';
	return
end
bb=0;
for i=1:32;
	bb=2*bb+b(i);
end
majorVersion=floor(bb/256/256/256);
minorVersion=bitand(15,floor(bb/256/256/16));
bugVersion=bitand(15,floor(bb/256/256));
release=bitand(255,floor(bb/256));
nonRelRev=bitand(255,floor(bb));
switch(release/16)
	case 2; x=sprintf('d%x',nonRelRev);
	case 4; x=sprintf('a%x',nonRelRev);
	case 6; x=sprintf('b%x',nonRelRev);
	case 8;
		if nonRelRev>0
			x=sprintf('f%x',nonRelRev);
		else
			x='';
		end
	otherwise
		error(sprintf('The Gestalt string ''%s'' does not yield a standard Apple version code.',gestaltString));
end
if bugVersion~=0
	versionString=sprintf('%d.%d.%d%s',majorVersion,minorVersion,bugVersion,x);
else
	versionString=sprintf('%d.%d%s',majorVersion,minorVersion,x);
end
