function WaitTicks(wait)
% WaitTicks(wait)
%
% Wait the requested number of system ticks. One tick is 1/60.15 seconds.
%
% TIMING ADVICE: the first time you access any MEX function or M file,
% Matlab takes several hundred milliseconds to load it from disk.
% Allocating a variable takes time too. Usually you'll want to omit
% those delays from your timing measurements by making sure all the
% functions you use are loaded and that all the variables you use are
% allocated, before you start timing. MEX files stay loaded until you
% flush the MEX files (e.g. by changing directory or calling CLEAR
% MEX). M files and variables stay in memory until you clear them.
%
% 
% Windows: ________________________________________________________________
%
% WaitTicks is implemented using WaitSecs.m, and each tick is interpreted 
% as 1/60.15 seconds to be consistent with the Mac version.
% (In Windows, a system tick is usually 1 millisecond, but with a precision that
% varies from system to system. We're ignoring the system ticks here.)
%
% OS X: ___________________________________________________________________
%
% "Ticks" functions are deprecated, use "Secs" functions instead;  To
% measure time, use GetSecs insted of GetTicks. To delay, use WaitSecs
% instead of WaitTicks.  
%
% "Secs" functions are provided for OS 9, Windows, and OS X and are more precise
% and more accurate than "Ticks" functions.  
%
% WaitTicks behaves exactly as on OS 9, relying on the Psychtoolbox mex
% function GetTicks to read the system tick count.
% _________________________________________________________________________
%
% See also: GetTicks, GetSecs, WaitSecs.
%
% 6/14/95 dhb  Added to help.
% 1/29/97 dhb  More comments.
% 3/15/97 dgp  Expanded comments.
% 3/15/99 xmz  Put in conditional for Windows.
% 2/12/04 awi  Added the OS X case.  Fixed the Windows case so that it
%              waits an integer number of ticks.  Simplified platform test
%              expressions.  Added test for unfamiliar platform.  Added OS X 
%              section to comments.  





if IsWindows
    wait=ceil(wait);
    WaitSecs(wait/60.15);
elseif IsOSX
    ticks=GetTicks+wait;
	while GetTicks<ticks
         
    end
else
    error('Platform unrecognized by Psychtoolbox');
end  %else    
    
    
    
        
