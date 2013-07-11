function MovieDemo
% MovieDemo
% 
% Show a movie, easy as pie.
%
% The comments below review the difference between the OS 9 version of
% the toolbox and PTB-3.
% 
% New and different on OS Psychtoolbox: 
% 
%  The OS 9 Psychtoolbox is based on QuickTime and PTB-3 is
%  based on OpenGL.  Differences between these underlying graphics libraries
%  have unavoidably resulted in differences between the OS 9 and PTB-3
%  versions of Screen; Some subfunctions of OS 9 Screen are replaced by
%  similar, but not identical functions in PTB-3.
%
%  In OS 9:                   In PTB-3
%   
%   WaitBlanking              Flip        
%   PutImage                  MakeTexture
%   CopyWindow                DrawTexture
% 
%  WaitVBL vs. Flip
%
%   In OS 9 an onscreen window consisted of only one buffer: the front
%   buffer.  Any pixel drawn to an onscreen window would automatically
%   appear on the display as soon as the monitor scan passed  that pixel,
%   the maximum delay between drawing a pixel and its appearnce on the
%   monitor being one frame period.  To  prevent vertical shearing of
%   animated displays, Psyhctoolbox scripts called WaitBlanking before
%   updating the display.  WaitBlanking synchonrized drawing into vide RAM
%   wih the monitor scan by delaying the onset of drawing until the scan
%   reached the bottom right corener of the display, at which point a
%   vertical blank would occur while the CRT beam retraced to the upper
%   left corner of the display.
% 
%   In contrast, PTB-3 onscreen windows by default have not one but two
%   buffers:  the front buffer and the back buffer.  In double-buffered
%   mode, all Screen drawing commands are issued to the back buffer, where
%   what is drawn is hidden from view.  The contents of the back buffer
%   remain hidden from view until The Screen('Flip', ...) command is
%   issued. Flip waits  until the vertical retrace and then interchages the
%   contents of the front and back buffers, bringining into view on the
%   display what had previosly been hidden in the back buffer. 
%
%  OpenOffscrenWindow + PutImage + CopyWindow vs. MakeTexture + DrawTexture
%
%   In OS 9, PutImage copied a MATLAB matrix into an offsceen window.  The
%   offscreen Window could then be quickly copied to the display during an
%   animation.  Thus  PutImage was one part of three for quickly
%   displaying MATLAB matrices during an animation:
%
%       Steps to quickly display an image matrix in OS 9:
%         1. Create an offscreen window by calling OpenOffscreenWindow
%         2. Copy the matrix to an offsceen window using PutImage
%         3. Copy the offscreen window to an onscreen window using
%         CopyWindow
%
%   Instead of using Offscreen Windows to quickly display MATLAB Matrices 
%   the OSX Psychtoolbox uses Textures.
%
%       Steps to quickly display an image matrix in PTB-3:  
%         1. Create a texture and copy into it a matrix by using MakeTexture
%         2. Copy the texture quickly to the onscreen window using
%         DrawTexure
%
%   What is the difference betwen offscreen windows and textures? The
%   difference is that Drawing commands such as FillOval may be  issued to
%   offscreen windows but not to Textures; You may only fill textures by
%   using matrices. 
%
%   How do we draw ovals during an animation if we can not store them
%   within offscren or textures? Unlike QuickDraw, OpenGL is fast enought
%   to render all drawing commands directly to the display during an
%   animation loop. Except for MATLAB matrices stored in textures, there
%   should be no need to prerender and buffer what is displayed during
%   animation.  
%
%  Animation loops: OS 9 VS PTB-3
%
%   To understand the differences between OS 9 and PTB-3 Screen its helpful
%   to compare animation loops:
%
%   OS 9: 
%
%       %create and fill offscreen windows here... 
%       for i=1:numberOfMovieFrames
%           Screen(window, 'WaitBlanking');
%           Screen('CopyWindow', offscreenWindow(i), window);
%       end
%
%   PTB-3: 
%
%       %generate textures here...
%       for i=1:numberOfMovieFrames
%           Screen('DrawTexture', texture(i), window);
%           Screen(window, 'Flip');
%       end
%
%   
% See also: DriftDemo, PsychDemos.
%
%
% HISTORY
%
% 7/3/04   awi  Wrote it.  Based on Denis Pelli's MovieDemo for OS 9.
% 7/19/04  awi  restored Priority, and time tests, added option to plot results. 
% 7/24/04  awi  Cosmetic
% 9/8/04   awi  Added Try/Catch.
% 4/23/05  mk   Added Priority(0) to Catch-Section.
% 11/19/06 dhb  Remove OSX from name.


try
	fprintf('MovieDemo.m is short. Read it to see how easy it \nis to show movies.\n\n');
	
	% If there are multiple displays guess that one without the menu bar is the
	% best choice.  Dislay 0 has the menu bar.  
	screens=Screen('Screens');
	screenNumber=max(screens);
	
	% Open a window.  Note the new argument to OpenWindow with value 2,
	% specifying the number of buffers to the onscreen window.  
	[window,windowRect]=Screen(screenNumber,'OpenWindow',0,[],[],2);
	
	
	% Give the display a moment to recover from the change of display mode when
	% opening a window. It takes some monitors and LCD scan converters a few seconds to resync.
	WaitSecs(2);
	
	
	% Show the movie, first forwards, then backwards.  We don't need offscreen
	% windows as in OS 9.  FillOval is now fast enough to draw ovals during the
	% animation.
	n=300;
	rect=[0 0 n n];
	rect2=AlignRect(rect,windowRect,RectRight,RectBottom);
	
	black=BlackIndex(window);
	white=WhiteIndex(window);
	Screen(window, 'FillRect', white);
	Screen(window,'TextColor', black);
	%Screen(window, 'TextFont', 'Courier');
	%Screen(window,'TextSize',24);
	for(i=1:2)
        Screen(window, 'FillRect', white);
        Screen(window,'DrawText','Showing movie at priority 0 ...',10,30,black);
        Screen(window,'Flip');
	end
	WaitSecs(2)
	t=0;
	HideCursor;
	for i=[1:n/2, n/2:-1:1]
		r=[0 0 2 2]*(i);
        Screen(window,'FillRect', white);
		Screen(window,'FillOval',black,r);
        t=t+1;
        s(t)=Screen(window, 'Flip');
	end
	ShowCursor;
	s1Diffs=diff(s);
	
	Screen(window, 'FillRect', white);
	Screen(window,'DrawText','Showing movie at priority 9 ...',10,30,black);
	Screen(window,'Flip');
	WaitSecs(2);
	priorityLevel=MaxPriority(window);
	t=0;
	HideCursor;
	Priority(priorityLevel);
	%this is strange, raising the priority seems to cause the first call to
	%flip to return immediatly.  Needs investigation.  
	Screen(window, 'Flip');
	for i=[1:n/2, n/2:-1:1]
		r=[0 0 2 2]*(i);
        Screen(window,'FillRect', white);
		Screen(window,'FillOval',black,r);
        t=t+1;
        s(t)=Screen(window, 'Flip');
	end
	Priority(0);
	ShowCursor;
	s2Diffs=diff(s);
	
	
	% In OS 9 we checked for missed frames by comparing the actual movie play
	% duration to the expected duration as predicted by the number of frames
	% shown and the frame rate  reported by Screen('FrameRate'). However,
	% FrameRate on PTB-3 returns the nominal frame rate reported by the video
	% card driver.  That is imprecise.  So instead, we use the median frame
	% period, which will work except for in case of egregious synching failure.
	frameRate=1/median(s2Diffs);
	frames1=sum(s1Diffs)*frameRate-length(s1Diffs);
	frames2=sum(s2Diffs)*frameRate-length(s2Diffs);
	
	% Close the on- and off-screen windows
	Screen('CloseAll');
	
	if round(frames1)==0 && round(frames2)==0
		s=sprintf('Success! The movie was shown twice. Both showings were frame-accurate.');
	else
		s=sprintf('The movie was shown twice, running over by %.0f frames in the first and %.0f frames in the second showing.',frames1,frames2);
	end
	s=sprintf('%s The first showing may have been jerky, due to interruptions. The second showing used Priority, at priority %g, to minimize interruptions, to make the movie run smoothly. See Priority and Rush.',s,priorityLevel);
	fprintf('%s\n',WrapString(s));
	
	doPlotInput=input('Plot frame intervals [y/n]? ','s');
	if streq(upper(doPlotInput), 'Y') || strcmpi(doPlotInput, 'YES')
        plot(s2Diffs);
        xlabel('frame number');
        ylabel('interval (seconds)');
	end %if
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    ShowCursor;
    Priority(0);
    psychrethrow(psychlasterror);
end % try..catch
