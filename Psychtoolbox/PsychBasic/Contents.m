% 
%   Psychtoolbox:PsychBasic
%   Psychtoolbox:PsychBeta:PsychBasic
%  
%   help Psychtoolbox % For an overview, triple-click me & hit enter.
%   help PsychDemos   % For demos, triple-click me & hit enter.
%
%     Beeper               - Play a nice beep tone of selectable duration, frequency and volume.
%     CharAvail            - Is a keypress available for GetChar?       
%     DisableKeysForKbCheck - Tell KbCheck and KbWait to ignore specific keys.
%     DoNothing            - Does nothing. Used to time Matlab's overhead.
%     DrawFormattedText    - Drawing of formatted text into windows.
%     FlushEvents          - Flush any unprocessed events. 
%     FontInfo             - Return a struct array describing installed fonts.
%     Gestalt              - Query system configuration on OS 9 and OS X. 
%     GetBusTicks          - Number of system bus ticks since startup.
%     GetBusTicksTick      - Duration of one tick of the GetBusTicks clock.
%     GetChar              - Wait for keyboard character and return it.
%     GetPID               - Get the process ID of the MATLAB process.
%     GetMouse             - Get mouse position. 
%     GetMouseWheel        - Get mouse wheel position delta on a wheel mouse.
%     GetSecs              - Time since startup with high precision. 
%     GetSecsTick          - Duration of one tick of the GetSecs clock.
%     GetTicks             - Number of 60.15 Hz ticks since startup. 
%     GetTicksTick         - Duration of one tick of the GetTicks clock.
%     HideCursor           - Hide cursor.
%     IOPort               - A I/O driver for access to serial ports.
%     KbCheck              - Get instantaneous keyboard state.
%     KbEventAvail         - Return number of pending keyboard events in ringbuffer.
%     KbEventFlush         - Remove all pending keyboard events in ringbuffer.
%     KbEventGet           - Get oldest pending keyboard event in ringbuffer.
%     KbKeysAction         - Return an incremented or decremented value, depending on keys pressed.
%     KbName               - Convert keycode to key name and vice versa.
%     KbPressWait          - Wait for key press, make sure no keys pressed before.
%     KbQueueCreate        - Create keyboard queue.
%     KbQueueRelease       - Destroy keyboard queue.
%     KbQueueFlush         - Empty keyboard queue.
%     KbQueueStart         - Start recording of key presses into queue.
%     KbQueueStop          - Stop recording of key presses into queue.
%     KbQueueCheck         - Check keyboard queue for key presses/releases.
%     KbReleaseWait        - Wait until all keys on keyboard are released.
%     KbStrokeWait         - Wait for single, isolated key stroke.
%     KbTriggerWait        - Wait for trigger keys on keyboard.
%     KbWait               - Wait until at least one key is pressed and return its time.
%     ListenChar           - Start GetChar queue.
%     LoadPsychHID         - Helper function for loading PsychHID on MS-Windows.
%     MachAbsoluteTimeClockFrequency - Mach Kernel time measurement.  
%     PredictVisualOnsetForTime - Predict stimulus onset for given Screen('Flip') 'when' timespec.
%     psychassert          - Drop in replacement for Matlabs assert().
%     psychlasterror       - Drop in replacement for Matlabs lasterror().
%     psychrethrow         - Drop in replacement for Matlabs rethrow().
%     PsychCV              - Miscellaneous C routines for computer vision and related stuff.
%     PsychDrawSprites2D   - Fast drawing of many 2D sprite textures.
%     PsychKinect          - Psychtoolbox driver for the Microsoft XBOX-360 Kinect.
%     PsychtoolboxDate     - Current version date, e.g. '1 August 1998'
%     PsychtoolboxVersion  - Current version number, e.g. 2.32
%     PsychWatchDog        - Watchdog mechanism and error handler for Psychtoolbox.
%     PsychTweak           - Tweak Psychtoolbox low-level operating parameters.
%     RemapMouse           - Map mouse position to stimulus position.
%     RestrictKeysForKbCheck - Restrict operation of KbCheck et al. to a subset of keys on the keyboard.
%     Screen               - Control the video display. ** Type "Screen" for a list. ** 
%     SetMouse             - Set mouse position.
%     ShowCursor           - Show the cursor, and set cursor type.
%     Snd                  - Play sounds.
%     VideoRefreshFromMeasurement - Alternative calibration procedure to find exact video refresh interval.
%     WaitSecs             - Wait specified time.
%     WaitTicks            - Wait specified number of 60.15 Hz ticks.

% from OS 9
%
%   Psychtoolbox:PsychBasic.
%  
%   help Psychtoolbox % For an overview, triple-click me & hit enter.
%   help PsychDemos   % For demos, triple-click me & hit enter.
%  
%     Bytes               - How much memory is free? (MEX)
%     CharAvail           - Is a keypress available for GetChar? (Uses EventAvail.mex)
%     Debugger            - Enter low-level debugger. (MEX)
%     DoNothing           - Does nothing. Used to time Matlab's overhead. (MEX)
%     EventAvail          - Check for events: mouse, keyboard, etc. (MEX)
%     FileShare           - Control filesharing. (Uses FS.mex)
%     FlushEvents         - Flush any unprocessed events. (MEX)
%     FrameRate           - Quick accurate old measurement of frame rate, in Hz.
%     GetChar             - Wait for keyboard character and return it (with time & modifiers).
%     GetClicks           - Wait for mouse click(s); get location and number. (MEX)
%     GetMouse            - Get mouse position. (MEX)
%     GetSecs             - Time since startup (20 us precision, or better). (MEX)
%     GetTicks            - Number of 60.15 Hz ticks since startup. (MEX)
%     HideCursor          - Hide cursor. (MEX)
%     KbCheck             - Get instantaneous keyboard state. (MEX, fast)
%     KbName              - Convert keycode to key name.
%     KbWait              - Wait for key press and return its time. (MEX, fast)
%     LoadClut            - Loads the CLUT, supports all DAC sizes and pixelSizes.
%     MaxPriority         - The maximum priority compatible with a list of functions. 
%     PatchTrap           - Disable Mac OS routines. Dangerous! (MEX)
%     PrepareScreen       - Called by Screen when first window is opened on a screen.
%     Priority            - Disable interrupts. Dangerous! (MEX)
%     PsychtoolboxDate    - Current version date, e.g. '1 August 1998'
%     PsychtoolboxVersion - Current version number, e.g. 2.32
%     RestoreScreen       - Called by Screen when last window is closed (among all screens) or Screen.mex is flushed.
%     Rush                - Execute code quickly, minimizing interrupts. (MEX)
%     Screen              - Fifty display functions. ** Type "screen" for a list. ** (MEX)
%     ScreenSaver         - Control screen savers, eg AfterDark, Sleeper. (MEX)
%     SetMouse            - Set mouse position. (MEX)
%     ShowCursor          - Show the cursor, and set cursor type. (MEX)
%     Showtime            - Create and show QuickTime movies. (MEX)
%     Snd                 - Play sounds. (MEX)
%     WaitSecs            - Wait specified time. (MEX)
%     WaitTicks           - Wait specified number of 60.15 Hz ticks.

