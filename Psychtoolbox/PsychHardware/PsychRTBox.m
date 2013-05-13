function varargout = PsychRTBox(varargin)
% Driver for the USTC reaction time button box (RTBox) by Xiangrui Li et al.
% varargout = PsychRTBox(cmd, varargin);
%
% This driver allows to control most functions of the USTC RTBox response
% button box. In theory this driver should support boxes up to Box/Firmware
% version 5. In practice it has only been tested by the Psychtoolbox
% developer up to firmware and box version 1.3, therefore you may encounter
% bugs for later versions. The box itself comes bundled with an alternative
% driver called "RTBox" which is maintained by the developers of the box
% hardware itself, in case this driver doesn't work with your box or
% firmware.
%
% The RTBox is a USB device which provides 4 response buttons (pushbuttons)
% for subject responses and can report any button press- or release by the
% subject. Additionally it has an input for reporting of external
% electronic trigger signals and a photo-diode input for reporting of
% visual stimulus onset. The box uses a built-in high-resolution clock to
% timestamp all button- or trigger events, independent of the host
% computers clock in order to make it more reliable for response time
% measurements than most other response devices. It also buffers all events
% internally, so experiment scripts can read back events when it is most
% convenient. Timestamps can be either reported in Psychtoolbox standard
% GetSecs timebase for direct comparison with timestamps from GetSecs,
% WaitSecs, KbCheck et al., Screen('Flip') and PsychPortAudio, etc. This
% simplifies reaction time calculations. Timestamps can also be reported in
% the timebase of the boxe, e.g., time of a button press relative to the
% photo-diode light trigger signal or electronic trigger signal, if this is
% more convenient for a given experiment setup.
%
% Current versions of the RTBox have additional functionality, e.g., more
% digital trigger inputs, some sound trigger, and TTL trigger outputs.
%
% See http://lobes.usc.edu/RTbox for up to date product information and
% additional driver software.
%
% Please note that while the device documentation claims that the external
% electronic pulse port is able to receive TTL trigger signals, we couldn't
% verify that this is the case with our test sample of version 1 of the
% RTbox hardware. While the pulse input port responded to some pulses sent
% by one TTL compatible device, it failed to detect the majority of signals
% from TTL most other test devices.
%
% This indicates that the pulse port may not be fully TTL compliant and
% will need additional tinkering for your setup. However, results with
% later versions of the hardware may be different than our experience with
% our test sample.
% 
%
% The following subcommands are currently suppported:
% ===================================================
%
%
% handle = PsychRTBox('Open' [, deviceID] [, skipSync=0]);
% -- Try to open a connected RTBox, return a device handle 'handle' to it
% on success. The handle can be used in all further subcommands to refer to
% the box. By default, all USB ports (or rather USB-Serial ports) are scanned
% for a connected RTBox and the driver will connect to the first box found.
% Alternatively you can specify which box to use via the optional
% 'deviceID' namestring. This can be either the name of a box, or the name
% of the USB-Serial port to which the box is connected. This way you can avoid
% scanning of all ports and disambiguate in case multiple boxes are
% connected to your computer.
%
% Btw., if you only make use of one single RTBox, you don't need to specify
% the 'handle' parameter to all following subfunctions. Instead you can
% specify that parameter as [] or omit it and the driver will use the only
% open box connected.
%
% The optional parameter 'skipSync', if set to 1, will prevent the open
% routine from performing an initial clock synchronization. By default, it
% will perform an initial clock synchronization.
%
% After opening the box, you may want to invoke this method:
%
%
% clockRatio = PsychRTBox('ClockRatio' [, handle] [, durationSecs]);
% -- Perform a clock drift calibration between the computers GetSecs host
% clock and the internal clock of the box 'handle'. Restrict calibration to
% a maximum of 'durationSecs' (default 60 seconds if omitted). Return the
% computed 'clockRatio' and use it for all further operations.
%
% Due to manufacturing imperfections and environmental factors, no two
% clocks ever run at exactly the same speed. Therefore the computer clock
% and box clock will slowly "drift out of sync" under normal conditions,
% rendering retrieved event timestamps inaccurate over the course of a long
% experiment session. This calibration routine will exercise the clocks and
% compute the clock drift due to this speed difference, then use the
% computed drift (= clockRatio) to correct all reported timestamps for this
% drift, thereby providing the accuracy needed for reaction time studies.
%
% The clockRatio value tells, how many seconds of GetSecs time elapse when
% the box clock measures 1 second elapsed time. Ideally this value would be
% 1, ie. both clocks run at the same speed. A more realistic value would be,
% e.g., 1.000009 -- The computer clock goes 9 microseconds faster than the
% box clock, so the drift will accumulate an error of 9 microseconds for
% each elapsed second of your study.
%
% As every calibration, this routine involves some measurement/calibration
% error and is therefore not perfect, so even after a successfull
% 'ClockRatio' calibration, timestamps reported during your experiment will
% accumulate some error during the course of long experiment sessions.
%
% There are multiple ways to handle this:
%
% a) Use a long calibration time in this function for accurate results, and
% a reasonably short experiment duration.
%
% b) Repeat this procedure after every large block of trials, ie., every
% couple of minutes, e.g., while the subject is allowed to take a break in
% a long experiment session.
%
% c) Use the PsychRTBox('SyncClocks') function after each short block of
% trials, or even after each trial, for the highest accuracy.
%
% d) Don't care for clock drift throughout the experiment session, just
% collect the event timestamps in box clock format (see the 3rd return
% argument of PsychRTBox('GetSecs') or the returned timing array of
% PsychRTBox('BoxSecs');) and store them in some array. Remap all timestamps
% into the computers GetSecs time at the end of your session via
% PsychRTBox('BoxsecsToGetsecs'). This requires a bit more
% discipline from you in programming and organizing your data, but it
% provides the most accurate timestamps.
%
%
% [syncResult, clockRatio] = PsychRTBox('SyncClocks' [, handle]);
% -- Synchronize or resynchronize the clocks of the host computer and the
% box. Return result in 'syncResult' and the current clockRatio in
% 'clockRatio'. This routine is automatically carried out during invocation
% of PsychRTBox('ClockRatio'); but you can repeat the sync procedure
% anytime between trials via this subfunction for extra accuracy at the
% expense of about 0.5 - 1 second additional time for each invocation. You
% would typically execute this function at the start of each large block of
% trials, or before start of each trial if you are really picky about
% super-exact timing. The syncResult contains three values:
%
% syncResult(1) = Host time (GetSecs time) at time of clock sync.
%
% syncResult(2) = Box time at time of clock sync.
%
% syncResult(3) = Confidence interval for the accuracy of the sync. This
% value (in seconds) provides a reliable upper bound for the possible error
% introduced in all reported timestamps from the box. The real error may be
% significantly smaller, this is just an upper bound that you can check.
% Typical results on a well working system should be in the sub-millisecond
% range, e.g., 0.0003 seconds or 0.3 msecs. Typical results on a rather
% noisy system would be around 0.001 second or 1 msec. Results worse than 2
% msecs indicate some problem with your system setup that should be fixed
% before executing any experiment study which involves reaction time
% measurements. By default, the sync procedure will abort with an error if
% it can't calibrate to an accuracy with a maximum error of 1.3 msecs within
% a duration of 0.5 seconds. You can change these default constraints with
% a call to PsychRTBox('SyncConstraints').
%
%
% [oldmaxDurationSecs, oldgoodEnoughSecs, oldrequiredSecs, oldsyncMethod] = PsychRTBox('SyncConstraints'[, maxDurationSecs][, goodEnoughSecs][, requiredSecs][, syncMethod]);
% -- Change the constraints to apply during calls to PsychRTBox('SyncClocks');
% Optionally return old settings.
%
% 'maxDurationSecs' limits any call to 'SyncClocks' to a duration of at
% most the given number of seconds. Calibration aborts after at most that
% time, even if unsuccessfull - in that case with an error message. By
% default, the duration is limited to 0.5 seconds.
% 'goodEnoughSecs' Calibration will finish before 'maxDurationSecs' have
% elapsed, if the result is more accurate than an error of at most
% 'goodEnoughSecs'. By default, this is set to zero seconds, i.e.,
% calibration will always take 'maxDurationSecs'.
% 'requiredSecs' - The calibration will only use samples with an
% uncertainty of at most 'requiredSecs'. If not even a single sample of the
% required precision can be acquired within 'maxDurationSecs', the call
% will fail with an error, indicating that your system setup doesn't
% provide the required timing precision for your demands. By default, the
% minimum required precision is 0.0013 seconds, ie., it will tolerate an
% error of at most 1.3 msecs.
% 'syncMethod' - Select the synchronization method to use. Either, method 0
% (prewrite sync), or method 1 (postwrite sync), or method 2 (average). If
% you want to know the difference between the methods, please consult the
% source code of this file and read the code for the subroutine 'function
% syncClocks'. All three methods are robust and accurate within the
% returned confidence window, usually better than 1 msec. So far, method 1
% seems to get the best results on our test setups, so this is the default
% for the current driver release. However, we are still evaluating if
% method 0 would be a tiny bit better and worth switching the default to
% that.
%
%
% oldverbose = PsychRTBox('Verbosity' [, handle], verbosity);
% -- Set level of verbosity for driver: 0 = Shut up. 1 = Report errors
% only. 2 = Report warnings as well. 3 = Report additional status info. 4 =
% Be very verbose about what is going on. The default setting is 3 --
% Report moderate status output.
%
%
% devinfo = PsychRTBox('BoxInfo' [, handle] [, newdevinfo]);
% -- Return a struct 'devinfo' with all information about the current
% status and parameter settings for RTBox 'handle'. Optionally set a new
% struct with updated parameters via 'newdevinfo'. This function is mostly
% useful for debugging and benchmarking the driver itself. Most information
% contained in 'devinfo' will be useless for your purpose.
%
%
% PsychRTBox('Close', handle);
% -- Close connection to specific box 'handle'. Release all associated
% ressources.
%
%
% PsychRTBox('CloseAll');
% -- Close connections to all attached RTBox devices. Reset the PsychRTBox
% driver completely. You'll usually use this function at the end of your
% experiment script to clean up.
%
%
% oldeventspec = PsychRTBox('Enable' [,handle][, eventspec]);
% -- Enable specified type of event 'eventspec' on box 'handle'. This
% allows to enable detection and reporting of a specific type of event. By
% default, only reporting of push-button press is enabled, as this is the
% most common use of a response box.
%
% The following names are valid for the name string 'eventspec':
% 'press' = Report push-button press. This is the default setting.
% 'release' = Report push-button release.
% 'pulse' = Report electronic trigger events on external input port.
% 'light' = Report reception of light flashes by photo-diode on light port.
% 'tr' = Report reception of scanner trigger "TR" (TTL input from pin 7 of DB-9 port).
% 'all' = Enable all events.
%
% If called without the 'eventspec' parameter, the function will return the
% names of all currently enabled events.
%
%
% oldeventspec = PsychRTBox('Disable' [,handle][, eventspec]);
% -- Disable specified type of event 'eventspec' on box 'handle'. This
% allows to disable detection and reporting of a specific type of event. By
% default, only reporting of push-button press is enabled, as this is the
% most common use of a response box.
%
% See 'Enable' call for help on parameters.
%
%
% Once you have setup and calibrated the box and selected the type of
% events to detect and report, you will want to actually retrieve
% information about events. For this you use these commands:
%
%
% PsychRTBox('Start' [, handle] [, dontwaitforstart=0]);
% -- Start event detection and reporting by the box. The box will start
% detecting button and trigger events from here on and record them in the
% event buffer.
% 
% You will usually call this at the beginning of a response period. By
% default, the box has reporting already enabled after 'Open'ing it.
%
% The optional 'dontwaitforstart' parameter, if set to 1, will ask the
% 'Start' function to return control as soon as possible, ie., without
% waiting for confirmation of the box that event reporting has actually
% started. By default, the routine waits for an acknowledgement from the
% box, which can take 16 - 30 msecs in some cases.
%
%
% PsychRTBox('Stop' [, handle]);
% -- Stop event detection and reporting by the box. The box will ignore
% detecting button and trigger events from here on and no longerrecord them
% in the event buffer.
% 
% You will usually call this at the end of a response period.
%
%
% PsychRTBox('Clear' [, handle] [, syncClocks=0] [, dontRestart=0]);
% -- Stop event detection and reporting by the box, clear all recorded
% events so far, then restart reporting if it was active before calling
% this function.
% 
% Instead of calling 'Start' and 'Stop' to mark the start and end of a
% response period in a trial you can also simply use this function at the
% beginning of a trial (or its response period) to discard any stale data
% from a previous trial (or non-response interval).
%
% You can prevent an automatic restart of event reporting by setting the
% optional flag 'dontRestart' to a value of 1.
%
% You can ask the box to resynchronize its clock to the host computer clock
% by setting the optional flag 'syncClocks' to a value of 1. This is the
% same as calling PsychRTBox('SyncClocks').
%
%
% [time, event, boxtime] = PsychRTBox('GetSecs' [, handle] [, interTimeout=0.1] [, maxTimeout=interTimeout] [, maxItems=inf]);
% -- Retrieve recorded events from the box 'handle'.
%
% By default, as many events are returned as are available within the
% test interval, but you can select a specific number of wanted events
% by setting the optional parameter 'maxItems'. If there aren't any pending
% events from the box, by default the driver waits for up to 0.1 seconds
% for events to arrive. You can change this 'interTimeout' interval via the
% positive (non-zero) 'interTimeout' parameter. The function will return if
% no new events show up within 'interTimeout' seconds. If something shows
% up, the deadline for return is extended by 'interTimeout' seconds. You
% can set an absolute upper limit to the response interval via the
% 'maxTimeout' parameter. That defaults to 'interTimeout' if omitted.
% Please note that after an event is detected by the box, up to 16-32 msecs
% can elapse until the event is received by the computer, so you may not
% want to set these timeout values too small!
%
% The function will return an array of timestamps in 'time', and an array
% of corresponding names of the events in 'event'. E.g., event(1) will
% report the identity of the first detected event, e.g., '1' if button 1
% was pressed, whereas time(1) will tell you when the event happened, ie.,
% when button 1 was pressed. 'time' is expressed in host clock time, aka
% GetSecs() time. If no events are pending since last invocation of this
% function, empty vectors will be returned.
%
% Additionally, the vector 'boxtime' contains the same timestamp, but
% expressed in box clock time. See below for a use of that.
%
% By default, the following names are possible for 'event's:
%
% '1' = 1st button pressed, '1up' = 1st button released.
% '2' = 2nd button pressed, '2up' = 2nd button released.
% '3' = 3rd button pressed, '3up' = 3rd button released.
% '4' = 4th button pressed, '4up' = 4th button released.
% 'pulse' = electronic pulse received on electronic pulse input port.
% 'light' = Light pulse received by photo-diode connected to light input port.
% 'tr' = Scanner trigger "TR" (TTL input from pin 7 of DB-9 port) received.
% 'serial' = PsychRTBox('Trigger') Softwaretrigger signal received on USB-Serial port.
%
% Note: 'tr' is only supported on boxes with Firmware version 3.0 or later.
%
% However, you can assign arbitrary names to the buttons and events if you
% don't like this nomenclature via the PsychRTBox('ButtonNames') command.
%
% The reported timestamps are expressed in host clock time, ie., in the
% same units as the timestamps returned by GetSecs, Screen('Flip'),
% PsychPortAudio, KbCheck, KbWait, etc., so you can directly calculate
% reaction times to auditory stimuli, visual stimuli and other events.
%
% See the help for PsychRTBox('SyncClocks') and PsychRTBox('ClockRatio')
% for the accuracy of these timestamps and tips for obtaining optimal
% accuracy.
%
% Additionally the event times are also returned in 'boxtime', but this
% time expressed in box time -- the time of the box internal clock.
% 
%
% There are multiple variants of this query command with the same optional
% input arguments, but different return arguments. All of these return
% timestamps in box time without remapping to GetSecs time by calling:
%
% [boxtime, event] = PsychRTBox('BoxSecs' ...);
% -- Timestamps are in raw box clock time, everything else is the same as
% in PsychRTBox('GetSecs' ...).
%
% If you have the 'boxtime' timestamps from one of the previous functions
% around, you can map them later to GetSecs time with very high precision
% at the end of your experiment session via:
%
% [GetSecs, Stddev] = PsychRTBox('BoxsecsToGetsecs' [, handle], boxTimes);
% -- Perform a post-hoc mapping of a vector of raw box timestamps
% 'boxTimes' into a vector of host clock 'GetSecs' timestamps. Return some
% error measure in 'Stddev' as well, if available.
%
% This method can be used to convert event timestamps expressed in the box
% clocks timebase into timestamps in Psychtoolbox GetSecs host clock
% timebase. It has the advantage of providing the highest possible accuracy
% in mapping, because it computes an optimal mapping function for this
% purpose, which is based on all the timing information collected
% throughout a whole experiment session. The disadvantage is that it will
% only provide meaningful results if you call it at the end of your
% experiment session, so you'll need to manage all your collected
% timestamps in a format that is suitable as input to this function.
%
%
% Timestamps can also be returned relative to a specific trigger event: You
% specify which event acts as a trigger. Then all timestamps of all events
% are expressed relative to the time of that trigger event, i.e., as
% deltas. Any event can be the trigger. Format of all arguments is
% as in PsychRTBox('BoxSecs' ...);
%
% E.g., PsychRTBox('serial', ...); Returns timestamps relative to the first
% occurence of a electronic input port trigger signal since the last query.
% PsychRTBox('light', ...); Returns timestamps relative to photo-diode
% light pulse. PsychRTBox('1'); returns relative to press of 1st button,
% etc. etc.
%
%
% sendTime = PsychRTBox('SerialTrigger' [, handle]);
% -- Send a software generated trigger to the box via the serial port
% connection. This will register as a event of type 'serial' and you can
% retrieve timestamps relative to the first trigger within a response
% period via the PsychRTBox('serial', ...); command.
%
%
% sendTime = PsychRTBox('EngageLightTrigger [, handle]);
% sendTime = PsychRTBox('EngagePulseTrigger [, handle]);
% sendTime = PsychRTBox('EngageTRTrigger [, handle]);
% 
% -- Engage trigger input on the box for reception of a one-shot trigger
% signal. This function will return immediately after submitting the
% request to the box. It may take up to 5 msecs worst-case until the
% trigger input is really enabled. If you want to wait for the trigger to
% be really enabled, call, e.g., PsychRTBox('Enable', handle, 'lighton'); instead,
% as that function will wait until the trigger is really active.
%
% Trigger events are special: If a trigger has been received, the
% box auto-disables the trigger input, preventing reception of any
% further trigger events, until the trigger gets reenabled. The trigger gets
% reenabled on many occasions if it has been enabled once via the
% PsychRTBox('Enable', ...); command, e.g., at each call to
% PsychRTBox('Start'); or PsychRTBox('Clear'). If you want to enable the
% trigger on-the-fly, then this function is your friend.
%
% The reason why light trigger auto-disables itself is because a typical
% CRT display monitor would generate such trigger signals at the rate of
% video refresh, once your stimulus is displayed, e.g., at a rate of 100
% Hz. Usually you only want to know one defined timestamp of initial
% stimulus onset, therefore the box prevents reception of all but the
% first light trigger.
%
% Similar reasoning applies to Pulse and TR triggers.
%
%
% oldNames = PsychRTBox('ButtonNames' [, handle] [, newNames]);
% -- Query or assign labels for the four response box buttons other than
% the default names.
%
% This function allows to assign arbitrary names to the four buttons on the
% box. These names are reported when querying for button presses and
% releases. By default, oldNames = PychRTBox('ButtonNames') would return
% the cell array with the four following names: '1', '2', '3', '4'. These
% are the names reported for button presses. Button releases would report
% the names with an 'up' appended, ie., '1up', '2up', '3up', '4up'. You can
% assign arbitrary new names by passing a cell array with four namestrings,
% e.g., PsychRTBox('ButtonNames', [], {'7', 'whats', 'hick', 'screw'})
% would assign the names '7', 'whats', 'hick' and 'screw' for button press
% events, and '7up', 'whatsup', 'hickup' and 'screwup' for release events
% of the corresponding buttons.
%
% Please note that the assignment of names to buttons must be unique, ie.
% assigning the same name to multiple buttons is not allowed.
%
%
% oldIntervals = PsychRTBox('DebounceInterval' [, handle] [, debounceSecs]);
% -- Query current button debounce intervals (in 4-element vector
% 'oldIntervals', one value for each button), and optionally set new
% debounce interval in seonds via the optional argument 'debounceSecs'.
% 'debounceSecs' can be a scalar, in which case the same setting is applied
% to all buttons, or a 4-element row vector, e.g., [0.1, 0.1, 0.1, 0.1] to
% set an individual interval for each of the four buttons.
%
% The built-in debouncer prevents multiple button responses (button press
% or release actions) from getting recorded/reported within some
% 'debounceSecs' debounce interval. After a button has changed state, only
% the type (press or release), identity (which button) and timestamp (when)
% of the first state change is reported. Any other state change within
% 'debounceSecs' seconds of time after that first change will be ignored.
% After that time has elapsed, further state changes are reported again. By
% default, this dead "debounce" interval is set to 0.050 seconds, ie., 50
% msecs. Button bouncing happens if a subject presses or releases a button
% very rapidly or vigorously. If such quick multiple events or bounces are
% not ignored, they will create multiple apparent button responses which
% are a hazzle to deal with in experiment scripts and data analysis.
%
% If you find multiple responses generated for only one apparent button
% press during piloting, you may want to set a bigger debounce interval
% with this function.
%
% Please note that debouncing doesn't apply to the PsychRTBox('ButtonDown')
% function.
%
% Please also note that there is another hardware debouncer with a duration
% of 0.3 msecs on RTBox versions with firmware versions 1.3 and older, which
% can't be disabled, so even if you'd set a zero interval here, you'd still
% get a minimum 0.3 msecs debounce period from the hardware itself.
%
% Later versions of the firmware support the following
% PsychRTBox('HardwareDebounce') command to control the hardware debounce
% interval more fine-grained.
%
%
% [oldValue] = PsychRTBox('HardwareDebounce' [, handle] [, scanNum]);
% -- Set/get hardware debouncer setting. The hardware will treat a button
% event as valid only if the button state stays stable for at least
% 'scanNum' scanning iterations of the firmware. The scan interval is about
% 67 microseconds. The valid scanNum is from 1 through 255, with a default
% setting of 16 cycles for 1.072 msecs debounce interval.
%
% For software debouncing at the driver level, see PsychRTBox('DebounceInterval')
% above. 
% 
%
% buttonState = PsychRTBox('ButtonDown' [, handle] [, whichButtons]);
% -- This reports the current button state of all response buttons of box
% 'handle', or a subset of response buttons if specified by the optional
% 'whichButtons' argument, e.g., whichButton = {'1', '4'} to only test
% buttons 1 and 4. 'buttonState' is a vector which contains a 1 for each
% pressed button, and a zero for each released button.
%
% This query is as instantaneous and "live" as possible. The reported state
% is not subject to button debouncing, but the measured "raw state".
% Usually you will want to use the PsychRTBox('GetSecs' ...) functions ans
% similar functions to query timestamped button state. They are typically
% as fast as this method and they provide timestamps of when the state was
% queried, whereas this function doesn't give you information about how
% "fresh" or recent the query is. However for simple button queries outside
% the response interval, e.g., while the box is PsychRTBox('Stop')'ped with
% no need for timestamps, this may be an option.
%
% Due to the design of the USB bus, the query may be outdated wrt. to the
% real state by up to 16 - 21 msecs, depending on operating system and
% driver configuration.
%
%
% buttonState = PsychRTBox('WaitButtonDown' [, handle] [, whichButtons]);
% -- Wait until at least one of the specified buttons in 'whichButtons' is
% pressed down. If 'whichButtons' is omitted, all buttons are tested.
%
%
% PsychRTBox('WaitButtonUp' [, handle] [, whichButtons]);
% -- Wait until all of the specified buttons in 'whichButtons' are
% released. If 'whichButtons' is omitted, all buttons are tested.
%
%
% [timeSent, confidence] = PsychRTBox('TTL' [, handle] [, eventCode=1]);
% - Send TTL to DB-25 port (pin 8 is bit 0). The second input is event code
% (default 1 if omitted), 4-bit (0~15) for box versions < 5, and 8-bit
% (0~255) for later versions. It can also be equivalent binary string, such
% as '0011'.
%
% The optional return arguments are the 'timeSent' when the TTL update was
% performed, and an upper bound on the uncertainty 'confidence' of
% 'timeSent'.
%
% The width (duration) of the TTL pulse is controlled by the
% PsychRTBox('TTLWidth') command.
%
% This function is only supported for v3.0 RTBoxes and later, the ones with
% EEG event code support.
%
% 
% [oldValue] = PsychRTBox('TTLWidth' [, handle][, widthSecs]);
% - Set/get TTL pulse width in seconds. The default width is 0.97e-3, ie.
% 97 microseconds when the device is opened. The actual width may have some
% small variation. The supported width ranges from 0.14e-3 to 35e-3 secs. A
% infinite width 'inf' is also supported. Infinite width means the TTL will
% stay until it is changed by the next PsychRTBox('TTL') command, such as
% PsychRTBox('TTL',0).
%
% This function is only supported for v3.0 RTBoxes and later, the ones with
% EEG event code port support.
%
% In Version <5.0, the TTL width at DB-25 pins 17~24 is controlled by a
% potentiometer inside the box. In Version >= 5, the width is also
% controlled by 'TTLWidth' command. 
% 
%
% [oldValue] = RTBox('TTLResting' [, handle][, newLevel]);
% - Set/get TTL polarity for DB-25 pins 1~8. The default is 0, meaning the
% TTL resting is low. If you set newLevel to nonzero, the resting TTL will
% be high level. If you need different polarity for different pins, let us
% know. This function is only supported with firmware version 3.1 and
% later.
% 
% In Version 5.0 and later, newLevel has second value, which is the
% polarity for pins 17~24.
% 

% TODO:
%
% - Debouncing for PsychRTBox('Buttondown') as well, or leave it "raw"?
%

% History:
% 08/01/2008 Initial implementation based on RTBox.m from Xiangrui Li (MK).
% 01/29/2009 "Close to beta" release. First checkin to SVN (MK).
% 01/30/2009 Improved syncClocks algorithm, option to spec a specific box
%            by port in the open call (MK).
% 02/08/2009 Huge redesign of API and internal routines. Now we use an
%            internal queue (MK).
% 02/14/2009 Refinements and cleanup (MK).
% 02/15/2009 More refinements and rework of post-hoc timestamp remapping (MK).
% 06/07/2009 Check for ambiguous assignment of buttonnames to avoid errors.
%            Bug found by Vinzenz Schoenfelder (MK).
% 06/14/2009 Remove special case code for Octave. No longer needed (MK).
% 12/14/2009 Update for RTBox'es with firmware v1.4 and later:
%
%            * Store box firmware version numerically for easier comparison.
%
%            * v >= 1.4 sends '?' acknowledge for live button query before
%              the byte that reports buttons state --> Handle this.
%
%            * Detect clock frequency of RTBox and adapt our bytes2secs() mapping
%              accordingly for V1.4 and later.
%
%            * 'lightoff' event is dead, now its called 'tr' instead for TR
%              scanner trigger reception on new V3.0 firmware boxes.
%
%            * Switch v1.3 and later boxes back to 'x' simple E-Prime et
%              al. compatible mode at PsychRTBox('close') time.
%
%            * Support new 'hardwaredebounce' command on v1.4+
%
%            * Support new 'ttlwidth' for v3.0, 'ttlresting' for v3.1+
%
%            * Support 4-bit TTL out port for v3.0+
%
%            * Add fast-engage commands for pulse and TR triggers as well.
%
% 10/20/2011 Switch fast calibration from robustfit to polyfit() - Now we
%            always use polyfit() in all cases to avoid need for Matlab
%            statistics toolbox. The PostHoc routine always did this, but
%            the online routine didn't. (MK)
%
% 11/26/2011 Updates for Firmware versions up to V5 (MK):
%
%            * TTL out is now 8 bit capable and 'ttlwidth' and 'ttlresting'
%              needed to be updated as well.
%            * Event enable handling has changed.
%            * Various other stuff.
%
%            None of this is tested due to lack of hardware/firmware.
%
% 01/06/2012 Bugfix for Firmware versions >= 4.1. Did not receive events
%            due to wrong acknowledgement handling. (MK)
%

% Global variables: Need to be persistent across driver invocation and
% shared with internal subfunctions:
global rtbox_info;
global rtbox_global;

% Start of driver code -- Entry point:

    % First time invocation? Perform init of device arrays and global
    % settings:
    if isempty(rtbox_info)
        % Setup device info struct array, as well as per-device default
        % settings:
        % CAUTION: Same settings are reassigned in the openRTBox()
        % subfunction each time PsychRTBox('Open') is called! The settings
        % made there override the settings made here!!!
        rtbox_info=struct('events',{{'1' '2' '3' '4' '1up' '2up' '3up' '4up' 'pulse' 'light' 'tr' 'serial'}},...
                              'enabled',[], 'ID','','handle',-1,'portname',[],'sync',[],'version',[],'clkRatio',1,'verbosity',3, ...
                              'busyUntil', 0, 'boxScanning', 0, 'ackTokens', [], 'buttons', [0 0 0 0; 0 0 0 0; 0 0 0 0], ...
                              'syncSamples', [], 'recQueue', [], 'boxClockTickIntervalSecs', 1/115200);

        % Setup event codes:
        rtbox_global.eventcodes=[49:2:55 50:2:56 97 48 57 89]; % code for 12 events
        
        % List of supported subcommands:
        rtbox_global.cmds={'close', 'closeall', 'clear', 'stop', 'start', 'test', 'buttondown', 'buttonnames', 'enable', 'disable', 'clockratio', 'syncclocks', ...
              'box2getsecs', 'boxinfo', 'getcurrentboxtime','verbosity','syncconstraints', 'boxsecstogetsecs', 'serialtrigger', ...
              'debounceinterval', 'engagelighttrigger', 'waitbuttondown', 'waitbuttonup', 'ttlwidth', 'ttlresting', 'hardwaredebounce', ...
              'ttl', 'engagepulsetrigger', 'engagetrtrigger' };
          
        % Names of events that can be enabled/disabled for reporting:
        rtbox_global.events4enable={'press' 'release' 'pulse' 'light' 'tr' 'all'};
        
        % Low-level protocol codes corresponding to the events:
        rtbox_global.enableCode='DUPOFA'; % char to enable above events, lower case to disable
        
        % Preload some functions of PTB we'll need:
        eval('GetSecs;WaitSecs(0.001);');
        
        % Selection of blocking strategy that the IOPort driver shall use
        % for blocking writes:
        if IsWin || IsOSX
            % A blocking wait gives very good results on OS/X. On Windows
            % it gives the same results as a polling wait (setting 2), so
            % we prefer blocking for lower cpu load at same quality on
            % Windows as well:
            rtbox_global.blocking = 1;
        else
            % On Linux, a polling wait is of advantage, so we use that:
            rtbox_global.blocking = 2;
        end

        % No devices open at first invocation:
        rtbox_global.nrOpen = 0;
        
        % Default settings for the syncClocks() function:
        % -----------------------------------------------
        
        % Use syncClocks() method 1 with postwrite timestamp by default:
        rtbox_global.syncmode = 1;

        % Maximum duration of a syncClocks calibration run is 0.5 seconds:
        rtbox_global.maxDuration = 0.5;

        % Desired 'minwin' calibration accuracy is 0.0 msecs: If we manage
        % to get better than that, we abort sampling. We also abort
        % sampling of the rtbox_global.maxDuration is reached:
        rtbox_global.optMinwinThreshold = 0.0;

        % Maximum allowable (ie. worst) acceptable minwin for a sample:
        % We default to 1.3 msecs, as a 1.2 msecs minwin is basically never
        % exceeded. It is unlikely that all samples within a syncClocks run
        % are worse than 1.3 msecs and the run would therefore fail.
        rtbox_global.maxMinwinThreshold = 0.0013;
        
        % Worst case delay after a command has been received by the box,
        % before it gets actually dequeued from the microprocessors serial
        % receive buffer and executed: 5 msecs is a very generous value to
        % be on the safe side:
        rtbox_global.maxbusy = 0.005;
    end

    if nargin < 1
        error('You must provide a command string to PsychRTBox!');
    end

    % Command dispatch:
    cmd = lower(varargin{1});
    if isempty(cmd)
        error('You must provide a non-empty command string to PsychRTBox!');
    end
    
    if strcmp(cmd, 'closeall') % Close all devices
        % Only close our devices, not other devices that may be opened
        % via IOPort but unrelated to us:
        for i=1:length(rtbox_info)
            s=rtbox_info(i).handle;
            if s>=0
                % Disable all scanning on box before close:
                stopBox(i);

                % Close connection:
                IOPort('Close', s);
                
                rtbox_info(i).handle = -1;
            end
        end

        rtbox_global.nrOpen = 0;       % Reset count of open devices to zero.
        clear rtbox_info;              % clear main device info struct array.
        clear rtbox_global;            % clear main global settings struct.
        
        return;
    end
    
    % Open the connection to device, do initial setup and sync:
    if strcmp(cmd, 'open')
        % Assign deviceID identifier of device to open, or the default
        % name 'Default' if none specified: User can also specify a serial
        % port device name for a device:
        if nargin < 2
            deviceID=[];
        else
            deviceID=varargin{2};
        end

        if isempty(deviceID)
            deviceID = 'Default';
        end

        if nargin < 3
            skipSync = [];
        else
            skipSync = varargin{3};
        end
        
        if isempty(skipSync)
            skipSync = 0;
        end
        
        % Open and initialize box:
        openRTBox(deviceID, rtbox_global.nrOpen+1);

        % Increment count of open boxes:
        rtbox_global.nrOpen = rtbox_global.nrOpen + 1;

        % Return as handle:
        varargout{1} = rtbox_global.nrOpen;

        if ~skipSync
            % Perform initial mandatory clock sync:
            syncClocks(rtbox_global.nrOpen);
        end
        
        % Perform initial button state query:
        buttonQuery(rtbox_global.nrOpen);
        
        % Start event scanning on box, with the above default enabled setting,
        % i.e., only button press 'D' reporting active:
        startBox(rtbox_global.nrOpen, 1);
        
        return;
    end

    if strcmp(cmd, 'syncconstraints')
        % Return current constraint settings:
        varargout{1} = rtbox_global.maxDuration;
        varargout{2} = rtbox_global.optMinwinThreshold;
        varargout{3} = rtbox_global.maxMinwinThreshold;
        varargout{4} = rtbox_global.syncmode;
        
        % Set constraints for syncClocks:
        if nargin > 1 && ~isempty(varargin{2})
            rtbox_global.maxDuration = varargin{2};
        end

        if nargin > 2 && ~isempty(varargin{3})
            rtbox_global.optMinwinThreshold = varargin{3};
        end

        if nargin > 3 && ~isempty(varargin{4})
            rtbox_global.maxMinwinThreshold = varargin{4};
        end
        
        if nargin > 4 && ~isempty(varargin{5})
            rtbox_global.syncmode = varargin{5};
        end

        return;
    end
    
    % Deal with variable number of inputs:
    if nargin - 1 > 0 
        nIn = nargin - 1;
    else
        nIn = 0;
    end
    
    if nIn > 1
        in2=varargin{3};
    end

    % Device handle provided?
    if nargin > 1
        % Yes: Store it in 'id':
        id = varargin{2};
    else
        % Nope.
        id = [];
    end

    % If no device handle - or empty default handle - provided, just
    % default to the first open RTBox device for convenience in
    % setups/scripts that only use one RTBox -- which is probably the
    % common case:
    if isempty(id)
        id = 1;
    end
    
    % Child protection:
    if ~isscalar(id) || id < 1 || id > length(rtbox_info) || id > rtbox_global.nrOpen
        error('Invalid device handle specified! Did you open the device already?');
    end

    if isempty(rtbox_info(id).handle)
        error('Invalid device handle specified! Maybe you closed this device already?');
    end

    % Build additional cell array of valid read commands:
    read=rtbox_info(id).events;
    read{end+1}='secs';    % Like GetSecs see below.
    read{end+1}='boxsecs'; % All events measured in absolute box time.
    read{end+1}='getsecs'; % All events measured in absolute GetSecs time.

    % Assign serial port handle:
    s = rtbox_info(id).handle;

    % Subcommand dispatch:
    switch cmd
        case 'verbosity'
            if nIn<2
                error('You must provide the new level of "verbosity" to assign!');
            end
            
            % Return old level of verbosity:
            varargout{1} = rtbox_info(id).verbosity;

            % Assign new level of verbosity to device:
            rtbox_info(id).verbosity = in2;
            
        case 'serialtrigger' % send serial trigger to device
            tWritten = sendTrigger(id);
            if nargout, varargout{1}=tWritten; end
            
        case 'engagelighttrigger'
            % Enable light on trigger quickly.
            tWritten = engageTrigger(id, 'O');
            if nargout, varargout{1}=tWritten; end
            
        case 'engagepulsetrigger'
            % Enable pulsetrigger quickly.
            tWritten = engageTrigger(id, 'P');
            if nargout, varargout{1}=tWritten; end
            
        case 'engagetrtrigger'
            % Enable TR-Trigger quickly.
            tWritten = engageTrigger(id, 'F');
            if nargout, varargout{1}=tWritten; end
            
        case 'ttl'
            % Send TTL 4 bit event to output port on supported hardware:
            if rtbox_info(id).version < 3
                RTBoxWarn('notSupported', in1, 3);
                return;
            end
            
            % Default event code is 1:
            if isempty(in2), in2 = 1; end

            % Can be a binary string:
            if ischar(in2), in2=bin2dec(in2); end
            
            % Range check:
            if rtbox_info(id).version < 5
                maxTTL = 15;
            else
                maxTTL = 255;
            end
            
            if (in2 < 0) || (in2 > maxTTL) || (in2~=round(in2))
                RTBoxError('invalidTTL');
            end
            
            % Decode to final trigger byte:
            if rtbox_info(id).version < 3.2
                in2 = dec2bin(in2,4);
                in2 = uint8(bin2dec(in2(4:-1:1))); % reverse bit order
            end

            if rtbox_info(id).version >= 5
                in2=[1 in2];
            end

            % Emit via blocking write:
            [tsend twin] = sendTTLPortEvent(id, in2);
            if nargout
                varargout={tsend twin};
            end

            if twin > 0.003
                fprintf('PsychRTBox: Warning! TTL trigger, send timestamp uncertainty %f msecs exceeds 3 msecs!\n', twin * 1000);
            end
            
        case 'debounceinterval'            
            % Return old debouncer settings for each button:
            varargout{1} = rtbox_info(id).buttons(3, :);

            if nIn<2
                return;
            end
            
            % Assign new settings:
            if isscalar(in2)
                % Single value: Apply to all buttons.
                rtbox_info(id).buttons(3, :) = [in2, in2, in2, in2];
            else
                % Multi value: Apply individually to each button:
                if size(in2,1)~=1 || size(in2, 2)~=4
                    error('Either set a single common debounce value for all buttons or pass a 4-element row vector with 4 settings for all 4 buttons!');
                end
                rtbox_info(id).buttons(3, :) = in2;                
            end
            
            % Reset debouncer:
            rtbox_info(id).buttons(2, :) = [0, 0, 0, 0];
            
        % Retrieve all pending events from the box, aka the serial port
        % receive buffers, parse them, filter/postprocess them, optionally
        % return mapped event timestamps in GetSecs timebase:
        case read % 12 triggers, plus 'secs' 'boxsecs' 'getsecs'
            cmdInd=strmatch(cmd,read,'exact'); % which command

            % Timestamp relative to trigger wanted?
            if cmdInd<13
                ind=[cmdInd<5 (cmdInd<9 && cmdInd>4) cmdInd==9:11];
                if ~rtbox_info(id).enabled(ind), RTboxError('triggerDisabled',rtbox_global.events4enable{ind}); end
                % minbytes=14; % at least 2 events
            end
            
            % Preinit return args to empty in case no event is detected:
            varargout={[] '' []};
            
            % 2nd argument is inter-response timeout: Return if there isn't
            % any data received for that amount of time. Each received new
            % item will extend that timeout by given amount:
            if nIn > 1 && ~isempty(varargin{3})
                intertimeout = varargin{3};
                if intertimeout <=0
                    error('Invalid interTimeout value specified. Must be significantly > 0 secs!');
                end
            else
                % Default is 0.1 secs aka 100 msecs:
                intertimeout = 0.1;
            end

            % 3rd argument is absolute timeout for all responses, the
            % absolute upper bound:
            if nIn > 2 && ~isempty(varargin{4})
                abstimeout = varargin{4};
                if abstimeout <=0
                    error('Invalid maxTimeout value specified. Must be significantly > 0 secs!');
                end
            else
                % Default is to set it to intertimeout:
                abstimeout = intertimeout;
            end
            
            % 4th argument is maximum number of responses to fetch at most:
            if nIn > 3 && ~isempty(varargin{5})
                maxItems = varargin{5};
            else
                % Default is to infinite, i.e., no limits: We get all we
                % can get within the setup timeout intervals:
                maxItems = inf;
            end
            
            % Retrieve events:
            [evid, timing] = getEvents(id, maxItems, maxItems, abstimeout, intertimeout);
            nevent = length(evid);
            
            % Anything retrieved?
            if nevent == 0
                return;
            end

            if nargout > 1 || cmdInd < 13
                % Map event id to human readable label string:
                for i=1:nevent % extract each event and time
                    ind=min(find(evid(i)==rtbox_global.eventcodes)); %#ok<MXFND> % which event
                    if isempty(ind)
                        RTboxWarn('invalidEvent',evid(i));
                        break; % not continue, rest must be messed up
                    end
                    event{i} = rtbox_info(id).events{ind}; %#ok event name
                end
            end

            if isempty(timing), return; end

            % Convert boxtiming and/or map it to host clock time:
            if cmdInd==15 || cmdInd==13
                % Convert into computer time: MK-Style

                % First return optional "raw" array with boxtimes:
                varargout{3} = timing;
                
                % Then remap to GetSecs host timebase:
                timing = box2GetSecsTime(id, timing);
                
            elseif cmdInd<13 % trigger: relative to trigger
                ind=strmatch(cmd,lower(event),'exact'); % trigger index
                if isempty(ind), RTboxWarn('noTrigger',cmd); return; end
                ind=ind(1); % use the 1st in case of more than 1 triggers
                trigT=timing(ind); % time of trigger event
                event(ind)=[]; timing(ind)=[]; % omit trigger and its time from output
                if isempty(event), return; end % if only trigger event, return empty
                timing=timing-trigT;   % relative to trigger time
            end

            varargout{1} = timing;

            if nargout > 1
                if length(event)==1, event=event{1}; end % if only 1 event, use string
                varargout{2} = event;
            end
            
        case 'boxinfo'
            % Return complete device info struct:
            varargout{1} = rtbox_info(id);
            
            % Optionally set a new one -- Only for debugging!!
            if nIn > 1
                rtbox_info(id) = in2;
            end
            
        case 'box2getsecs'
            % Map boxtime to GetSecs time with the recommended method:
            if nIn<2
                error('You must provide the boxtime to map!');
            end
            
            varargout{1} = box2GetSecsTime(id, varargin{3});
            
        case 'boxsecstogetsecs'
            % Map boxtime to GetSecs time post-hoc style:
            % We compute an optimal least-squares fit linear mapping
            % function of boxtime to hosttime, using all collected
            % syncClocks samples from the whole experiment session. Then we
            % remap all given input boxsecs samples to getsecs time via
            % lookup in that linear best-fit. This automatically corrects
            % for clock-drift and should provide the least possible error
            % on the mapping procedure, because it makes use of all
            % available sync information of a whole session:
            if nIn<2
                error('You must provide the boxtimes to map!');
            end
            
            [remapped, sd, clockratio] = box2GetSecsTimePostHoc(id, varargin{3}); 
            varargout = { remapped, sd , clockratio};
                                    
        case 'getcurrentboxtime'
            % Retrieve current time of box clock.
            
            % We do so by performing a syncClocks call, but retaining the
            % current .sync results, so we just (mis-)use the function for
            % our purpose:
            tmpsync = rtbox_info(id).sync;
            syncClocks(id); % clear buffer, sync clocks
            varargout{1} = rtbox_info(id).sync;
            rtbox_info(id).sync = tmpsync;

        case 'stop'
            % Stop event processing and reporting on box:
            % This will store all pending events in internal queue:
            stopBox(id);
            
        case 'start'
            % (Re-)Start event processing and reporting on box:

            if nIn >= 2 && ~isempty(varargin{3}) && varargin{3}
                % Asynchronous start, i.e. don't wait for acknowledge but
                % return asap:
                waitForStart = 0;
            else
                % Wait for acknowledge of start:
                waitForStart = 1;
            end
            
            startBox(id, waitForStart);
            
        case 'clear'
            % Clear all pending events on box, optionally perform clocksync:
            % By default, box is restarted after clear and no clocksync is
            % performed, but box is not restarted if usercode doesn't want
            % this or if it wasn't running before. Optionally clocksync is
            % executed:
            
            % Stop event processing on box, if active:
            boxActive = rtbox_info(id).boxScanning;
            if boxActive
                % This will store all pending events in internal queue:
                stopBox(id);
            end
            
            % Clear all buffers (serial buffers and event queue):
            purgeRTbox(id);
            
            % Optional syncClocks requested?
            if nIn >= 2 && ~isempty(varargin{3}) && varargin{3}
                % Perform clockSync:
                syncClocks(id);
                
                if nargout
                    varargout{1}=rtbox_info(id).sync;
                    varargout{2}=rtbox_info(id).clkRatio;
                end
            end
            
            % Restart box? We restart if it was running before and usercode
            % doesn't forbid a restart:
            if boxActive && (nIn < 3 || isempty(varargin{4}) || varargin{4} == 0)
                startBox(id, 1);
            end
            
        case 'syncclocks'
            % Synchronize host clock and box clock, i.e., establish mapping
            % between both:

            % Stop event processing on box, if active:
            boxActive = rtbox_info(id).boxScanning;
            if boxActive
                % This will store all pending events in internal queue:
                stopBox(id);
            end
            
            syncClocks(id); % clear buffer, sync clocks
            if nargout
                varargout{1}=rtbox_info(id).sync;
                varargout{2}=rtbox_info(id).clkRatio;
            end

            if boxActive
                % Restart box if it was running:
                startBox(id, 1);
            end
            
        case 'buttondown'
            % Perform query:
            b2 = buttonQuery(id);
            
            if nIn<2, in2=read(1:4); end % not specified which button
            in2=cellstr(in2); % convert it to cellstr if it isn't
            for i=1:length(in2)
                ind=strmatch(lower(in2{i}),read(1:4),'exact');
                if isempty(ind), RTboxError('invalidButtonName',in2{i}); end
                bState(i)=b2(ind); %#ok
            end

            varargout{1} = bState;
            
        case 'waitbuttondown'
            if nIn<2, in2=read(1:4); end % not specified which button
            in2=cellstr(in2); % convert it to cellstr if it isn't
            bState = zeros(1,length(in2));
            
            % Repeat queries until at least one of the target buttons is
            % down:
            while ~any(bState)
                b2 = buttonQuery(id);
                for i=1:length(in2)
                    ind=strmatch(lower(in2{i}),read(1:4),'exact');
                    if isempty(ind), RTboxError('invalidButtonName',in2{i}); end
                    bState(i)=b2(ind); %#ok
                end
            end
            
            varargout{1} = bState;

        case 'waitbuttonup'
            if nIn<2, in2=read(1:4); end % not specified which button
            in2=cellstr(in2); % convert it to cellstr if it isn't
            bState = ones(1,length(in2));
            
            % Repeat queries until all of the target buttons are up:
            while any(bState)
                b2 = buttonQuery(id);
                for i=1:length(in2)
                    ind=strmatch(lower(in2{i}),read(1:4),'exact');
                    if isempty(ind), RTboxError('invalidButtonName',in2{i}); end
                    bState(i)=b2(ind); %#ok
                end
            end
            
        case 'buttonnames' % set or query button names
            oldNames=rtbox_info(id).events(1:4);
            if nIn<2, varargout{1}=oldNames; return; end
            if length(in2)~=4 || ~iscellstr(in2), RTboxError('invalidButtonNames'); end
            
            for i=1:length(in2)
                in2{i} = lower(in2{i});
            end
            
            % Check for unambiguous assignment: Each buttonname must be
            % unique!
            for i=1:length(in2)
                if length(strmatch(lower(in2{i}), in2,'exact')) ~= 1
                    % Ambituous naming -> Same buttoname multiple times!
                    error('Same name assigned to multiple buttons! Sorry, names must be unique!');
                end
            end
            
            rtbox_info(id).events(1:4)=in2;
            for i=5:8
                rtbox_info(id).events(i)=cellstr([char(rtbox_info(id).events(i-4)) 'up']);
            end
            if nargout, varargout{1}=oldNames; end

            
        case 'ttlwidth'
            % One of the commands that change Firmware device settings.
            % Need to stop event processing and drain the event queue, as
            % the protocol doesn't implement proper acknowledge tokens:
            if rtbox_info(id).version < 3
                RTBoxWarn('notSupported', in1, 3);
                return;
            end
            
            % Stop event processing on box, if active:
            boxActive = rtbox_info(id).boxScanning;
            if boxActive
                % This will store all pending events in internal queue:
                stopBox(id);
            end

            % Perform active query for current firmware settings:
            b8=get8bytes(id);
            
            % TTL width unit in s, not very accurate:
            wUnit=0.139e-3; 
            
            % Return old/current setting as 1st argument:
            if rtbox_info(id).version > 4, b8(1) = 255 - b8(1); end
            varargout{1} = b8(1) * wUnit;

            % New settings provided?
            if nIn >= 2
                if isempty(in2)
                    % Default to 0.00097 seconds:
                    in2=0.00097;
                end

                % Infinity means: Disable TTL-width, hold setting until
                % manually changed:
                if in2 == inf
                    in2 = 0;
                end

                % Range check:
                if (in2 < wUnit * 0.9 || in2 > wUnit * 255 * 1.1) && (in2>0)
                    RTBoxWarn('invalidTTLwidth', wUnit);
                    return;
                end

                width = double(uint8(in2 / wUnit)) * wUnit;
                b8(1) = width / wUnit;
                if rtbox_info(id).version > 4, b8(1) = 255 - b8(1); end

                % Writeback new firmware settings:
                set8bytes(id, b8);
                
                if (in2 > 0) && (abs(width - in2) / in2 > 0.1)
                    RTBoxWarn('widthOffset', width);
                end
                
                if width==0
                    width=inf;
                end
                
                % Return new / current setting as 1st argument:
                varargout{1} = width;
            end
            
            if boxActive
                % Restart box if it was running:
                startBox(id, 1);
            end
            
        case 'ttlresting'
            % One of the commands that change Firmware device settings.
            % Need to stop event processing and drain the event queue, as
            % the protocol doesn't implement proper acknowledge tokens:
            if rtbox_info(id).version < 3.1
                RTBoxWarn('notSupported', in1, 3.1);
                return;
            end
            
            % Stop event processing on box, if active:
            boxActive = rtbox_info(id).boxScanning;
            if boxActive
                % This will store all pending events in internal queue:
                stopBox(id);
            end

            % Perform active query for current firmware settings:
            b8=get8bytes(id);
            
            % Return current setting as 1st return argument:
            varargout{1} = (b8(3) > 0);
            
            % New settings provided?
            if nIn >= 2
                if isempty(in2)
                    % Default to 0:
                    in2 = logical([0 1]);
                end
                                
                % Assign valid new setting:
                if rtbox_info(id).version < 5
                    b8(3) = in2(1) * 240; % '11110000'
                else
                    b8(3) = bitset(b8(3),1,in2(1));
                    b8(3) = bitset(b8(3),2,in2(2));
                end
                
                % Writeback new firmware settings:
                set8bytes(id, b8);
            end
                        
            if boxActive
                % Restart box if it was running:
                startBox(id, 1);
            end
            
            
        case 'hardwaredebounce'
            % One of the commands that change Firmware device settings.
            % Need to stop event processing and drain the event queue, as
            % the protocol doesn't implement proper acknowledge tokens:
            if rtbox_info(id).version < 1.4
                RTBoxWarn('notSupported', in1, 1.4);
                return;
            end
            
            % Stop event processing on box, if active:
            boxActive = rtbox_info(id).boxScanning;
            if boxActive
                % This will store all pending events in internal queue:
                stopBox(id);
            end

            % Perform active query for current firmware settings:
            b8=get8bytes(id);

            % Return current setting as 1st return argument:
            varargout{1} = b8(2);
            
            % New settings provided?
            if nIn >= 2
                if isempty(in2)
                    % Default to 16 firmware scancycles:
                    in2 = 16;
                end
                
                % Range check:
                if in2 < 1 || in2 > 255
                    % Invalid: Warn & Ignore:
                    RTBoxWarn('invalidScanNum');
                else
                    % Assign valid new setting:
                    b8(2) = uint8(in2);

                    % Writeback new firmware settings:
                    set8bytes(id, b8);
                end                
            end
            
            if boxActive
                % Restart box if it was running:
                startBox(id, 1);
            end
            
        case {'enable', 'disable'} % enable/disable event detection
            if nIn<2 % no event, return current state
                varargout{1}=rtbox_global.events4enable(rtbox_info(id).enabled);
                return;
            end

            % Stop event processing on box, if active:
            boxActive = rtbox_info(id).boxScanning;
            if boxActive
                % This will store all pending events in internal queue:
                stopBox(id);
            end
            
            isEnable=strcmp(cmd,'enable');
            in2=lower(cellstr(in2));
            for i=1:length(in2)
                ind=strmatch(in2{i},rtbox_global.events4enable,'exact');
                if isempty(ind), RTboxError('invalidEnable',rtbox_global.events4enable); end
                if ind==6, ind=1:5; end % all
                rtbox_info(id).enabled(ind)=isEnable; % update state
            end
            if nargout, varargout{1}=rtbox_global.events4enable(rtbox_info(id).enabled); end
            if ~any(rtbox_info(id).enabled), RTboxWarn('allDisabled',rtbox_info(id).ID); end
            
            if boxActive
                % Restart box if it was running:
                startBox(id, 1);
            end
            
        case 'clockratio' % measure clock ratio computer/box
            % Default to 60 seconds for clock ratio calibration, unless
            % specified otherwise:
            if nIn<2, in2=60; end;

            % Interval between trials
            interval = 1.5 * rtbox_global.maxDuration;
            
            % Number of trials to perform:
            ntrial=max(5,round(in2/interval));

            if rtbox_info(id).verbosity > 2
                fprintf('PsychRTBox: Measuring clock ratio on box "%s". Trials remaining:%4.f', rtbox_info(id).ID, ntrial);
            end
            
            % Stop event processing on box, if active:
            boxActive = rtbox_info(id).boxScanning;
            if boxActive
                % This will store all pending events in internal queue:
                stopBox(id);
            end
            
            % Switch to realtime priority if not already there:
            oldPriority=Priority;
            if oldPriority < MaxPriority('GetSecs')
                Priority(MaxPriority('GetSecs'));
            end
            
            % Perform ntrial calibration trials:
            tnow = GetSecs;
            for i=1:ntrial
                % Update rtbox_info.sync via a syncClocks() operation:
                syncClocks(id);
                
                % Store new syncClocks sample in array:
                t(i,:)=rtbox_info(id).sync; %#ok<AGROW>

                % Give status output:
                if rtbox_info(id).verbosity > 2
                    fprintf('\b\b\b\b%4.f',ntrial-i);
                end

                % Make sure that trials are at least 'interval' secs apart:
                tnow = WaitSecs('UntilTime', tnow + interval);
            end

            % Restore priority to old value:
            if Priority ~= oldPriority
                Priority(oldPriority);
            end
            
            % Delete syncClocks samples collected during this clockRatio
            % calibration from internal array, as they might bias later
            % box->host time remapping:
            rtbox_info(id).syncSamples = rtbox_info(id).syncSamples(1:end-ntrial, :);
            
            % Restart scanning on box if it was active before:
            if boxActive
                startBox(id, 1);
            end
            
            % Use always polyfit to fit a line (with least squares error)
            % to the samples of host clock and box clock measurements. We
            % use polyfit because it is part of default Matlab/Octave,
            % doesn't require statistics toolbox or fitting toolbox:
            [coef, st] = polyfit(t(:,2)-t(1,2),t(:,1)-t(1,1), 1);  % fit a line
            sd = st.normr * 1000; % std in ms

            % Assign new clock ratio for use by the timestamp mapping
            % routines later on:
            rtbox_info(id).clkRatio = coef(1);

            if rtbox_info(id).verbosity > 2
                fprintf('\n Clock ratio (computer/box): %.7f\n', rtbox_info(id).clkRatio);
                fprintf(' Sigma (ms): %.2g\n',sd);
            end
            
            if (sd > 0.5) && rtbox_info(id).verbosity > 1
                fprintf(2,' The ratio may not be reliable. Try longer time for ClockTest.\n');
            end

            % Optionally return new ratio:
            if nargout
                varargout{1} = rtbox_info(id).clkRatio;
            end

        case 'test' % for device test, no complain on error
            byte=IOPort('bytesAvailable',s);
            nevent=floor(byte/7);
            data=IOPort('read',s,1,nevent*7); % read data
            data=reshape(data,[7 nevent]);
            enableEvent(id, 'a'); % avoid interference to the following commands
            purgeRTbox(id); % clear buffer
            fprintf(' ID(%g): %s, %s\n',id,rtbox_info(id).ID,rtbox_info(id).version);
            fprintf(' Serial port: %s (handle = %g)\n',rtbox_info(id).portname, s);
            IOPort('write',s,'?'); % ask button state
            b2=IOPort('read',s,1,2);
            b2=dec2bin(b2(1)/16,4);
            b2=sprintf('%s ',b2(4:-1:1));
            fprintf(' Button down: %s\n',b2);
            fprintf(' Events enabled: %s\n',cell2str(rtbox_global.events4enable(find(rtbox_info(id).enabled)>0)));
            syncClocks(id); % sync clocks, restore detection
            fprintf(' ComputerClock/BoxClock: %.7f\n',rtbox_info(id).clkRatio);
            fprintf(' GetSecs-BoxSecs: %.4f\n',rtbox_info(id).sync(1));
            fprintf(' Number of events: %g\n',nevent);
            if nevent==0, return; end
            t0=bytes2secs(data(2:7,1), id); % first event time
            fprintf(' BoxSecs of event 1 (t0): %.3f s\n\n',t0);
            fprintf('  Data (7 bytes each)%19s  BoxSecs-t0\n','Events');
            for i=1:nevent  % disp each event and time
                ind=find(data(1,i)==rtbox_global.eventcodes); % which event
                if isempty(ind)
                    event=''; t=nan; % no error complain
                else
                    event=rtbox_info(id).events{ind}; % event ID
                    t=bytes2secs(data(2:7,i), id); % device time
                end
                fprintf('%4g',data(:,i));
                fprintf('%12s%8.3f\n',event,t-t0); % relative to first event
            end

        case 'close' % Close given device
            % Close port
            if ~isempty(rtbox_info(id))
                if s>=0
                    % Disable all scanning on box before close:
                    stopBox(id);

                    % Box with firmware version 1.3 or later?
                    if rtbox_info(id).version > 1.3
                        % Yes. This one supports "simple mode" for E-Prime
                        % compatibility as well, so switch it into that
                        % mode:
                        
                        % Enable all events
                        WaitSecs(0.1);
                        IOPort('Write', s, 'A');

                        % Switch to simple mode for E-Prime et al.:
                        WaitSecs(0.1);
                        IOPort('Write',s,'x');
                        
                        % Give some time to settle:
                        WaitSecs(0.1);
                    end
                    
                    % Close connection:
                    IOPort('Close', s);
                end

                % Mark corresponding struct slot as dead:
                rtbox_info(id).handle = [];

                % Decrease count of open devices:
                rtbox_global.nrOpen = rtbox_global.nrOpen - 1;
            end
            
        otherwise
            % Unknown command:
            RTboxError('unknownCmd',cmd,rtbox_global.cmds,rtbox_info(id),events);
            
    end % End of subfunction dispatch.

% End of drivers "main" routine:
end


% Helper functions:
% -----------------

% Map box timestamps to host clock (aka GetSecs) timestamps, based on clock
% sync results from syncClocks, and clock drift correction from clock ratio
% calibration:
function timing = box2GetSecsTime(id, timing)
    global rtbox_info;

    thostbase = rtbox_info(id).sync(1);
    tboxbase  = rtbox_info(id).sync(2);

    % Compute delta in box seconds from time when calibration was
    % done, wrt. box timer, then multiply with clockRatio
    % conversion factor to translate that delta into a delta wrt.
    % host clock, so timing will be the delta in host clock seconds
    % wrt. calibration time:
    timing = (timing - tboxbase) * rtbox_info(id).clkRatio;

    % Add basetime on host - this should give us correctly mapped
    % time in the GetSecs timebase of the computer:
    timing = timing + thostbase;
end

% Map box timestamps to host clock (aka GetSecs) timestamps, based on clock
% sync results from all syncClocks samples from a session:
function [timing, sd, clockratio] = box2GetSecsTimePostHoc(id, timing)
    global rtbox_info;
    global rtbox_global;

    % Stop box, if active, drain all queues. Please note that we don't
    % auto-restart the box, because this method is meant to be only called
    % at the end of an experiment session, so usercode must manually
    % restart if it wants to:
    stopBox(id);
    
    % Check if the latest syncClocks sample is older than 30 seconds. If
    % so, then we acquire a new final sample. We also resample if the last
    % sample is of too low accuracy, or if there are less than 2 samples in
    % total, as the fitting procedure needs at least 2 samples to work:
    while (size(rtbox_info(id).syncSamples, 1) < 2) || ...
          ((GetSecs - rtbox_info(id).syncSamples(end, 1)) > 30) || ...
          (rtbox_info(id).syncSamples(end, 3) > rtbox_global.maxMinwinThreshold)

      % Perform a syncClocks to get a fresh sample to finalize the sampleset:
      syncClocks(id);
    end

    % Extract samples for line fit:
    tbox  = rtbox_info(id).syncSamples(:, 2);
    thost = rtbox_info(id).syncSamples(:, 1);
    
    % MK: Change of strategy: We always use polyfit() instead of sometimes
    % regress() or robustfit() for the following reasons:
    %
    % 1. Strategy switching may introduce artifacts / non-predictable
    % behaviour into the results.
    %
    % 2. A simple polyfit() is better understood than the more advanced
    % robustfit() with its tons of tunable parameters and strategies, for
    % which i don't know which are appropriate for our data.
    %
    % 3. robustfit() and regress() require the Matlab statistics toolbox
    % and therefore won't work on a plain standard Matlab installation,
    % whereas polyfit() is part of default Matlab and Octave.
    %
    % The old code is left here for documentation purpose:
    %
    %     % More than 2 samples available?
    %     if 0 && length(tbox) > 2
    %         % Octave and older versions of Matlab don't have 'robustfit',
    %         % so we fall back to 'regress' if this function is lacking:
    %         if exist('robustfit') %#ok<EXIST>
    %             [coef st]=robustfit(tbox,thost);  % fit a line
    %             sd=st.robust_s; % stddev. in seconds.
    %         else
    %             coef =regress(thost, [ones(size(thost,1), 1), tbox ]);  % fit a line
    %             sd=0; % stddev. undefined with regress().
    %         end
    %         clockratio = coef(2);
    %         clockbias  = coef(1);
    %     else
    
    % Use polyfit, i.e., a least squares fit of a polynomial of either
    % degree 1 (purely linear), or degree 2 (quadratic). We choose a purely
    % linear fit if only 2 samples are available - the best we can do. If
    % more than 2 samples are available, we use a quadratic fit to account
    % for small non-linearities in the clock hardware over long intervals:
    [coef, st, mu] = polyfit(tbox, thost, min(length(tbox)-1, 2));

    clockratio = coef(1);
    clockbias  = coef(2); %#ok<NASGU>

    sd.norm = st.normr;

    %     end
    
    % Ok, got mapping equation getsecst = timing * clockratio + clockbias;
    % Apply it to our input timestamps:
    %    timing = timing * clockratio + clockbias;
    %
    %
    % Actually, we use polyval() to evaluate the equation. For linear case,
    % this is the same as above equation, but it conveniently also handles
    % the quadratic case and provides some margins on mapping error:
    if length(tbox) > 3
        [timing, sd.delta] = polyval(coef, timing, st, mu);
    else
        % sd.delta undefined for less than 4 elements in timing:
        sd.delta = repmat(Inf, size(timing));
        timing = polyval(coef, timing, st, mu);
    end

    if rtbox_info(id).verbosity > 3
        fprintf('PsychRTBox: In post-hoc box->host mapping: Coefficients of mapping polynom (in decreasing order) are: ');
        disp(coef);
        fprintf('\n\n');
    end
    
    % Ready.
end

% Clock sync routine: Synchronizes host clock (aka GetSecs time) to box
% internal clock via a sampling and calibration procedure:
function syncClocks(id)
    global rtbox_info;
    global rtbox_global;
    
    % Query level of verbosity to use:
    verbosity = rtbox_info(id).verbosity;
    
    % Perform multiple measurement trials per syncClocks run, take the best
    % one for final result. We use the "best" one because we have a good
    % criterion to find the best one.
    
    % We preallocate the sampling arrays for 250 samples at most. The
    % arrays will grow if 250 samples are not sufficient, at a small
    % performance penalty:
    ntrials = 250;
    
    % Any event reporting active?
    boxActive = rtbox_info(id).boxScanning;
    if boxActive
        % Disable all events on box. This will also clear all buffers:
        warning('PsychRTBox: syncClocks: Box was still scanning during invocation -- Driverbug?'); %#ok<WNTAG>
        stopBox(id);
    end

    % Switch to realtime priority if not already there:
    oldPriority=Priority;
    if oldPriority < MaxPriority('GetSecs')
        Priority(MaxPriority('GetSecs'));
    end
    
    % Get porthandle:
    s = rtbox_info(id).handle;
    t = zeros(3,ntrials);

    minwin = inf;
    tdeadline = GetSecs + rtbox_global.maxDuration;
    ic = 0;
    
    % Perform measurement trials until success criterion satisfied:
    % Either a sample with a maximum error 'minwin' less than desired
    % threshold, or maximum allowable calibration time reached:
    while (minwin > rtbox_global.optMinwinThreshold) && (GetSecs < tdeadline)
        % Wait some random fraction of a millisecond. This will desync us
        % from the USB duty cycle and increase the chance of getting a very
        % small time window between scheduling, execution and acknowledge
        % of the send operation:
        WaitSecs(rand / 1000);
        
        % Take pre-Write timestamp in tpre - Sync command not emitted
        % before that time. Write sync command, wait 'rtbox_global.blocking' for write
        % completion, store completion time in post-write timestamp tpost:
        [nw tpost, errmsg, tpre] = IOPort('Write', s, 'Y', rtbox_global.blocking);
        
        % We know that sync command emission has happened at some time
        % after tpre and before tpost. This by design of the USB
        % standard, host controllers and operating system USB stack. This
        % is the only thing we can take for granted wrt. timing, so the
        % "time window" between those two timestamps is our window of
        % uncertainty about the real host time when sync started. However,
        % on a well working system without massive system overload one can
        % be reasonably confident that the real emission of the sync
        % command happened no more than 1 msec before tpost. That is a
        % soft constraint however - useful for computing the final estimate
        % for hosttime, but nothing to be taken 100% for granted.
        
        % Write op successfull?
        if nw==1
            % Retrieve 7-Byte packet with timestamp from box:
            [b7 dummy errmsg] = IOPort('Read', s, 1, 7);
        else
            % Send op failed!
            fprintf('PsychRTBox: Warning! Sync token send operation to box failed! [%s]\n', errmsg);
            continue;
        end
        
        if length(b7)~=7 || b7(1)~=89
            % Receive op failed!
            fprintf('PsychRTBox: Warning! Corrupt data received from box! [%s]\n', errmsg);
            continue;
        end
        
        % Decode boxtime into seconds (uncorrected for clock-drift):
        %
        % From the computed raw box clock time in seconds, we subtract
        % another 0.00016 secs, aka 160 microseconds to compensate for the
        % average delay between reception of the sync token byte, and the
        % timestamping on the box. Delay is:
        % 1. Transmission delay from FTDI chip to Microprocessor over
        % serial link at 115200 baud with 10 bits (1 start + 8 data + 1
        % stop): 1000/115200*10 msecs = 0.087 msecs.
        % 2. Average scanning delay by firmware scanning loop: Interval is
        % 0.01 msecs, expected delay therefore 0.01/2 = 0.005 msecs.
        % 3. Time taken for firmware to take a snapshot of current box
        % clock time: 0.060 msecs.
        %
        % == 0.152 msecs. To account for minimal other delays, we round up
        % to 0.16 msecs aka 0.00016 secs and subtract this expected delay:
        tbox = bytes2secs(b7(2:7), id) - 0.00016;
        
        % Compute confidence interval for this sample:
        % For each measurement, the time window tpost - tpre defines a
        % worst case confidence interval for the "real" host system time
        % when the sync command was emitted.
        confidencewindow = tpost - tpre;
        
        % If the confidencewindow is greater than the maximum acceptable
        % window, then we reject this sample, else we accept it:
        if confidencewindow <= rtbox_global.maxMinwinThreshold
            % Within acceptable range. Accept this sample and check if it
            % is the best wrt. to window size so far:
            if confidencewindow < minwin
               % Best confidencewindow so far. Update minwin, as this is one
               % of the abortion critera:
               minwin = confidencewindow;
            end
            
            % Increase sample index to permanently accept this sample for
            % final set of competitors:
            ic = ic + 1;
            
            % Assign values:
            t(1,ic) = tpre;
            t(2,ic) = tpost;
            t(3,ic) = tbox;
        else
            % Inacceptably large error confidencewindow. Reject this sample:
            continue;
        end

        % Next sample pass:
    end

    % Done with sampling: We have up to 'ic' valid samples, unless minwin
    % is still == inf.
    
    % Restore priority to state pre syncClocks:
    if Priority ~= oldPriority
        Priority(oldPriority);
    end
    
    % Is box scanning supposed to be active?
    if boxActive
        % Restore event reporting:
        startBox(id, 1);
    end
    
    % At least one sample with acceptable precision acquired?
    if (minwin > rtbox_global.maxMinwinThreshold) || (ic < 1)
        % No, not even a single one! Set all results to "invalid"
        rtbox_info(id).sync=[inf, inf, inf];
        
        if verbosity > 1
            fprintf('PsychRTBox: Warning: On Box "%s", Clock sync failed due to confidence interval of best sample %f secs > allowable maximum %f secs.\n', rtbox_info(id).ID, minwin, rtbox_global.maxMinwinThreshold);
            fprintf('PsychRTBox: Warning: Likely your system is massively overloaded or misconfigured!\n');
            fprintf('PsychRTBox: Warning: See the help for PsychRTBox(''SyncConstraints'') on how to relax the constraints if you wish to do so.\n');
        end
        
        % Warn user:
        RTboxWarn('poorSync', minwin);
        
        % That's it:
        return;
    end

    % Ok, we have 'ic' > 0 samples with acceptable precision, according to
    % user specified constraints. Prune result array to valid samples 1 to ic:
    t = t(:, 1:ic);

    % We have three different methods for final sample selection. Each has
    % its favorable cases, but none of them is better by a large margin.
    % All can't guarantee accurate results on highly overloaded systems,
    % but only provide the best under given circumstances.
    %
    % Under low-load conditions, the new style methods may be slightly more
    % accurate under some conditions, but the difference is usually in the
    % sub-millisecond range, so either choice is ok. The options are mostly
    % here for benchmarking and stress testing of the driver...
    %
    % The default method (unless overriden by usercode or
    % userconfiguration) is method 0 -- New style with prewrite timestamps.
    % This because this method empirically seems to give best results on
    % low or medium load on our test set of machines and operating systems.
    switch rtbox_global.syncmode
        case 2
            % Old style method 2 - Middle of minwin confidence window:
            
            % Select sample with smallest confidence interval [tpre; tpost]:
            [mintdiff, idx] = min(t(2,:) - t(1,:));

            % Host time corresponds to midpoint of the confidence interval,
            % assuming a uniform distribution of likelyhood of the true write
            % time in the interval [tpre ; tpost]:
            hosttime = (t(1,idx) + t(2,idx)) / 2;
            
        case 1
            % New style method 1 - Postwrite timestamps:

            % Choose the most accurate sample from the set of candidates. This is
            % the sample with the smallest difference between the postwrite
            % timestamp and the associated box timestamp, ie., with the smallest
            % offset between postwrite host clock time and box clock time at
            % receive of sync command. The reasoning behind this goes like this:
            %
            % 1) The time offset between host clock and box clock is a constant -
            % at least within a significant multi-second time interval between
            % successive syncClocks calls (due to only small clock drift), but
            % certainly within a syncClocks run of a few hundred milliseconds
            % (error due to clock drift in this interval is negligible).
            %
            % 2) Computed box clock time t(3,:) is "close to perfect", as this
            % timestamp is taken by box microprocessor and firmware with a very
            % small and basically constant delay after sync token receive, ie.,
            % write completion. (Maximum theoretical error is smaller than 0.1 msecs).
            %
            % 3) The correct and optimal clock offset between host and box would be
            % tdiff = tsend - t(3,:) iff tsend would be host time at true write
            % completion.
            %
            % 4) The measured host time at write completion t(2,:) is always later
            % (and therefore numerically greater) than the true host time tsend at
            % write completion due to an unknown, random, greater than zero delay
            % tdelta, i.e., t(2,:) = tsend + tdelta, tdelta > 0. tdelta is the sum
            % of:
            %
            % a) The unknown delay of up to 1 msec between USB write-URB completion
            % by the USB host controller (which would be the real completion time
            % tsend) and detection of completion due to USB IOC (Interrupt-On-
            % Completion) due to invocation of the host controllers hardware
            % interrupt handler and host controller schedule scan and URB
            % retirement inside the interrupt handler.
            %
            % b) Random (and theoretically unbounded) scheduling delay / execution
            % delay between status update of the serial port data structures by the
            % interrupt handler and detection of write completion + timestamping by
            % the IOPort driver in polling mode, or scheduling delay between
            % wakeup-operation caused by the interrupt handler and start of
            % execution of the timestamping in the IOPort driver in blocking mode.
            %
            % The syncClocks error is therefore directly proportional to the size
            % of tdelta. Therefore:
            %
            % tdiff(:) = t(2,:) - t(3,:) by definition of clock offset host vs. box.
            % t(2,:) = tsend(:) + tdelta(:) by unknown scheduling/execution noise tdelta.
            %
            % It follows that by above definitions:
            %
            % tdiff(:) = tsend(:) - t(3,:) + tdelta(:);
            %
            % --> As we defined tsend(:) to be the unknown, but perfect and
            % noise-free, true send timestamp, and t(3,:) to be the perfect receive
            % timestamp by the box, it follows that by selecting the sample 'idx'
            % with the minimal tdiff(idx) from the set tdiff(:), we will select the
            % sample with the unknown, but minimal tdelta(idx). As tdelta accounts
            % for all the remaining calibration error, minimizing tdelta will
            % maximize the accuracy of the clock sync.
            %
            % ==> Select sample with minimum t(2,:) - t(3,:) as final best result:
            [mintdiff, idx] = min(t(2,:) - t(3,:));

            % mintdiff is our best estimate of clock offset host vs. box, and
            % t(:,idx) is the associated best sample. Unfortunately there isn't any
            % way to compute the exact residual calibration error tdelta(idx). The
            % only thing we know is that the error is bounded by the length of the
            % associated 'minwin' confidence interval of this sample, so we will
            % return 'minwin' as an upper bound on the calibration error. As
            % 'minwin' was used as a threshold in the sample loop for outlier
            % rejection, we can be certain that our estimate carries no greater
            % error than 'rtbox_global.maxMinwinThreshold'.

            % Extract all relevant values for the final sample:

            % Host time corresponds to tpost write timestamp, which should be as
            % close as possible to real host send timestamp:
            hosttime = t(2,idx);

        case 0
            % New style method 0 - Prewrite timestamps:

            % Choose the most accurate sample from the set of candidates.
            % This is the sample with the biggest difference between the
            % prewrite timestamp and the associated box timestamp, ie.,
            % with the biggest offset between prewrite host clock time and
            % box clock time at receive of sync command. The reasoning
            % behind this goes like this:
            %
            % 1) The time offset between host clock and box clock is a constant -
            % at least within a significant multi-second time interval between
            % successive syncClocks calls (due to only small clock drift), but
            % certainly within a syncClocks run of a few hundred milliseconds
            % (error due to clock drift in this interval is negligible).
            %
            % 2) Computed box clock time t(3,:) is "close to perfect", as this
            % timestamp is taken by box microprocessor and firmware with a very
            % small and basically constant delay after sync token receive, ie.,
            % write completion. (Maximum theoretical error is smaller than 0.1 msecs).
            %
            % 3) The correct and optimal clock offset between host and box would be
            % tdiff = tsend - t(3,:) where tsend would be unknown host time
            % at true write completion.
            %
            % 4) The measured host time before write submission t(1,:) is
            % always earlier
            % (and therefore numerically smaller) than the true host time tsend at
            % write completion due to an unknown, random, greater than zero delay
            % tdelta, i.e., t(1,:) = tsend - tdelta, tdelta > 0. tdelta is the sum
            % of:
            %
            % a) The unknown delay of up to 1 msec between USB write-URB completion
            % by the USB host controller (which would be the real completion time
            % tsend) and detection of completion due to USB IOC (Interrupt-On-
            % Completion) due to invocation of the host controllers hardware
            % interrupt handler and host controller schedule scan and URB
            % retirement inside the interrupt handler.
            %
            % b) Random (and theoretically unbounded) scheduling delay / execution
            % delay between status update of the serial port data structures by the
            % interrupt handler and detection of write completion + timestamping by
            % the IOPort driver in polling mode, or scheduling delay between
            % wakeup-operation caused by the interrupt handler and start of
            % execution of the timestamping in the IOPort driver in blocking mode.
            %
            % The syncClocks error is therefore directly proportional to the size
            % of tdelta. Therefore:
            %
            % tdiff(:) = t(1,:) - t(3,:) by definition of clock offset host vs. box.
            % t(1,:) = tsend(:) - tdelta(:) by unknown execution noise tdelta.
            %
            % It follows that by above definitions:
            %
            % tdiff(:) = tsend(:) - t(3,:) - tdelta(:);
            %
            % --> As we defined tsend(:) to be the unknown, but perfect and
            % noise-free, true send timestamp, and t(3,:) to be the perfect receive
            % timestamp by the box, it follows that by selecting the sample 'idx'
            % with the maximal tdiff(idx) from the set tdiff(:), we will select the
            % sample with the unknown, but minimal tdelta(idx). As tdelta accounts
            % for all the remaining calibration error, minimizing tdelta will
            % maximize the accuracy of the clock sync.
            %
            % ==> Select sample with maximum t(1,:) - t(3,:) as final best result:
            [maxtdiff, idx] = max(t(1,:) - t(3,:));

            % maxtdiff is our best estimate of clock offset host vs. box, and
            % t(:,idx) is the associated best sample. Unfortunately there isn't any
            % way to compute the exact residual calibration error tdelta(idx). The
            % only thing we know is that the error is bounded by the length of the
            % associated 'minwin' confidence interval of this sample, so we will
            % return 'minwin' as an upper bound on the calibration error. As
            % 'minwin' was used as a threshold in the sample loop for outlier
            % rejection, we can be certain that our estimate carries no greater
            % error than 'rtbox_global.maxMinwinThreshold'.

            % Extract all relevant values for the final sample:

            % Host time corresponds to tpre write timestamp, which should be as
            % close as possible to real host send timestamp:
            hosttime = t(1,idx);
        otherwise
            error('PsychRTBox: syncClocks: Unknown timestamping method provided. This is a driverbug!!');
    end
    
    % Box timers time taken "as is":
    boxtime  = t(3,idx);
    
    % Recalculate upper bound on worst case error 'minwin' from this best
    % samples tpost - tpre:
    minwin = t(2,idx) - t(1,idx);

    if verbosity > 3
        fprintf('PsychRTBox: ClockSync(%i): Box "%s": Got %i valid samples, maxconfidence interval = %f msecs, winner interval %f msecs.\n', rtbox_global.syncmode, rtbox_info(id).ID, ic, 1000 * rtbox_global.maxMinwinThreshold, 1000 * minwin);
        if rtbox_global.syncmode == 1
            fprintf('PsychRTBox: Confidence windows in interval [%f - %f] msecs. Range of clock offset variation: %f msecs.\n', 1000 * min(t(2,:)-t(1,:)), 1000 * max(t(2,:)-t(1,:)), 1000 * range(t(2,:) - t(3,:)));
        else
            fprintf('PsychRTBox: Confidence windows in interval [%f - %f] msecs. Range of clock offset variation: %f msecs.\n', 1000 * min(t(2,:)-t(1,:)), 1000 * max(t(2,:)-t(1,:)), 1000 * range(t(1,:) - t(3,:)));
        end
    end
    
    % Assign (host,box,confidence) sample to sync struct:
    rtbox_info(id).sync=[hosttime, boxtime, minwin];
    
    % Also store the sample in the syncSamples history:
    rtbox_info(id).syncSamples(end+1, :) = [hosttime, boxtime, minwin];
end

% Send enable/disable string to tell box to enable or disable detection and
% reporting of certain events. Wait for acknowledge from box that the
% requested event enable/disable was carried out properly. Certain
% characters encode certain types of events. A capital character enables
% that event type, a minor character disables that event type. The special
% letter 'A' or 'a' enables or disables all events. See the init routine
% for a mapping table of button/electronictrigger,lightTrigger,serialTrigger
% events to control characters.
%
% It is useful to disable event detection/reporting both to filter out
% unwanted events for a given experiment, but also during clock calibration
% and clockSync to prevent unwanted events from messing up the
% communication protocol between the box and the driver and from
% introducing unwanted processing delays during this time critical
% procedure.
function enableEvent(handle, str)
    global rtbox_info;
    global rtbox_global;

    if rtbox_info(handle).boxScanning
        error('PsychRTBox: enableEvent() called while box is scanning! Driverbug?'); %#ok<WNTAG>
    end
    
    s = rtbox_info(handle).handle;
    
    if isnumeric(str) && isscalar(str)
        % Firmware version >= 4.1 : Encodes enabled events by a single byte.
        
        % First send event setup command code "e", then the given enable
        % byte:
        str = [uint8('e') , uint8(str)];
    end
    
    for ie=1:length(str)
        % Wait until we can be sure that the box is ready to receive new
        % commands. The deadline is computed so that in the worst conceivable
        % case the box will be ready to receive at least 1 new command byte:
        WaitSecs('UntilTime', rtbox_info(handle).busyUntil);

        % Send control character for event enable/disable:
        [nw, tpost] = IOPort('Write', s, str(ie));

        % Command submission completed at time 'tpost'. Set a new busyUntil
        % time. After that time, the box should have stopped at the latest:
        rtbox_info(handle).busyUntil = tpost + rtbox_global.maxbusy;

	% On firmware versions < 4.1, we expect an acknowledge byte for each
        % sent byte. On later versions we only get an acknowledge for the
        % first 'e' character, but not the actual enable info byte afterwards:
	if (rtbox_info(handle).version < 4.1) || (ie == 1)
	    % Store this command code as one of the codes to be waited for
            % in parseQueue:
            rtbox_info(handle).ackTokens(end+1) = str(ie);
	end
    end

    % All enable tokens submitted.
    return;
end

% Retrieve firmware internal settings:
function b8 = get8bytes(handle)
    global rtbox_info;
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;
    
    % Write query command blocking:
    IOPort('Write', s, 's');
    
    % Wait blocking for 8 byte response:
    b8=IOPort('Read',s,1,8);
    
    return;
end

% Set/Update firmware internal settings:
function set8bytes(handle, b8)
    global rtbox_info;
    global rtbox_global;
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;

    % Write 'S'et command blocking:
    IOPort('Write',s,'S');

    % Wait blocking for 1 byte ack:
    IOPort('Read', s, 1, 1);
    
    % Write blocking new 8 byte settings:
    [nw tpost] = IOPort('Write', s, uint8(b8));

    % Command submission completed at time 'tpost'. Set a new busyUntil
    % time. After that time, the box should have completed command at the latest:
    rtbox_info(handle).busyUntil = tpost + rtbox_global.maxbusy;
    
    return;
end

function [tsend, twin] = sendTTLPortEvent(handle, ecode)
    global rtbox_info;
    global rtbox_global;
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;
    
    % Wait until we can be sure that the box is ready to receive new
    % commands. The deadline is computed so that in the worst conceivable
    % case the box will be ready to receive at least 1 new command byte:
    WaitSecs('UntilTime', rtbox_info(handle).busyUntil);

    % Submit event code to RTbox:
    [nw tpost err tpre] = IOPort('Write', s, ecode); % send

    % Compute tsend time approximation and confidence window twin:
    twin = tpost - tpre;
    tsend = tpre;
    
    % Command submission completed at time 'tpost'. Set a new busyUntil
    % time. After that time, the box should have completed command at the latest:
    rtbox_info(handle).busyUntil = tpost + rtbox_global.maxbusy;

    return;
end

% Send a 'Y' serial trigger command to box:
function tpost = sendTrigger(handle)
    global rtbox_info;
    global rtbox_global;
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;
    
    % Wait until we can be sure that the box is ready to receive new
    % commands. The deadline is computed so that in the worst conceivable
    % case the box will be ready to receive at least 1 new command byte:
    WaitSecs('UntilTime', rtbox_info(handle).busyUntil);
    
    % Box ready to receive our "serial trigger" command code
    % 'Y'. Send it blocking:
    [nw, tpost] = IOPort('Write', s, 'Y');
    
    % Command submission completed at time 'tpost'. Set a new busyUntil
    % time. After that time, the box should have stopped at the latest:
    rtbox_info(handle).busyUntil = tpost + rtbox_global.maxbusy;

    return;
end

% Send a "enable optical/pulse/TTL trigger" command evType to box:
function tpost = engageTrigger(handle, evType)
    global rtbox_info;
    global rtbox_global;
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;
    
    % Wait until we can be sure that the box is ready to receive new
    % commands. The deadline is computed so that in the worst conceivable
    % case the box will be ready to receive at least 1 new command byte:
    WaitSecs('UntilTime', rtbox_info(handle).busyUntil);
    
    % Box ready to receive our enable light trigger command code
    % evType. Send it blocking:
    [nw, tpost] = IOPort('Write', s, evType);
    
    % Command submission completed at time 'tpost'. Set a new busyUntil
    % time. After that time, the box should have stopped at the latest:
    rtbox_info(handle).busyUntil = tpost + rtbox_global.maxbusy;

    return;
end

% Query current state of box buttons -- the raw state, unaffected by
% debouncing or actual enable event settings, as close to now-time as
% possible:
function bState = buttonQuery(handle)
    global rtbox_info;
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;
    
    % Box in scanning mode? And button state change event reporting active?
    if (rtbox_info(handle).boxScanning) && all(rtbox_info(handle).enabled(1:2))
        % Yes. Box provides button press/release events, so its sufficient
        % to parse the receive queue.
        
        % Drain the serial receive buffer until either 10 seconds total
        % time elapsed or no new data for at least 20 msecs:
        parseQueue(handle, inf, 10, 0.020);
    else
        % Box is maybe scanning, but not reporting.
        
        % Is box scanning?
        isActive = rtbox_info(handle).boxScanning;
        if isActive
            % Yes. We can't perform a synchronous query with box in
            % scanning mode. Need to stopBox():
            stopBox(handle);
            
            % Now we startBox(), because this will trigger a callback into
            % us recursively, i.e. to buttonQuery(), but this time with
            % isActive == false, therefore the recursive call will go
            % through the else-clause and execute the query. -- No need for
            % us to do it here redundantly :-)
            startBox(handle, 1);
            
            % Done.
        else
            % Box is not scanning and reporting, all serial buffers are empty
            % and idle. We perform a synchronous query:
            IOPort('write',s,'?'); % ask button state: '4321'*16 63
            b2=IOPort('read',s, 1, 2); % ? returns 2 bytes
            if length(b2)~=2 || ...
               ((rtbox_info(handle).version  < 1.4) && (b2(2)~=63 || mod(b2(1),16)~=0)) || ...
               ((rtbox_info(handle).version >= 1.4) && (b2(1)~=63 || mod(b2(2),16)~=0))
                   % Failed! Corrupt response:
                   warning('PsychRTBox: Corrupt 2-byte response received in explicit buttonQuery()!'); %#ok<WNTAG>
            else
                % Uppermost 4 bits are button states:
                rtbox_info(handle).buttons(1, :) = [ 0, 0, 0, 0 ];
                % Set corresponding entries for pressed buttons:
                if rtbox_info(handle).version  < 1.4
                    % Version prior 1.4: Data first, then ?
                    rtbox_info(handle).buttons(1, logical(bitget(b2(1), 5:8))) = 1;
                else
                    % Version >= 1.4: ? first, then data:
                    rtbox_info(handle).buttons(1, logical(bitget(b2(2), 5:8))) = 1;
                end
            end
        end
    end
    
    bState = rtbox_info(handle).buttons(1, :);
    return;
end

% Start event scanning, detection and reporting on box:
function startBox(handle, WaitForStart)
    global rtbox_info;
    global rtbox_global;
    
    % Box scanning?
    if rtbox_info(handle).boxScanning
        % Box already started. Nothing to do:
        return;
    end
    
    % Emit a buttonstate query command, so we get a button state update:
    buttonQuery(handle);
    
    % Enable all events that are selected by usercode:
    customEnable = rtbox_info(handle).enabled(1:6);
    
    % Synchronize the button up/down events. Either both on or both off.
    % We do this so our software based button state live reporting in
    % parseQueue() works correctly - It needs press and release to be
    % reported:
    if any(customEnable(1:2))
        customEnable(1:2) = [1 1];
    else
        customEnable(1:2) = [0 0];
    end
    
    % The box always starts and stops with all event processing disabled,
    % transitions are only possible in startBox() and stopBox(), and we can
    % only enter startBox() after a previous call to stopBox() or at
    % initial startup. Therefore we know that all event processing is
    % disabled when we reach this point.
    %
    % enableEvent() will enable the user defined events non-blocking, ie.
    % without waiting for acknowledge from the box:
    
    if rtbox_info(handle).version < 4.1
        % Old style enable codes:
        enableEvent(handle, rtbox_global.enableCode(find(customEnable>0)) ); %#ok<FNDSB>
    else
        % Firmware version >= 4.1 : Encodes enabled events by a single byte.
        % Translate our customEnable into a single byte:
        byteEnable = 2.^(0:5) * customEnable';
        enableEvent(handle, byteEnable);
    end
    
    % Box event scanning and reporting is active:
    rtbox_info(handle).boxScanning = 1;
    
    % Should we wait until we get acknowledge from box?
    if WaitForStart
        % Yes. Stay in parseQueue until all event enable codes
        % acknowledged:
        parseQueue(handle, -2, inf, 10);
    end
    
    return;
end

% Stop event detection by box and all data transmission. Drain serial
% receive queue completely after box is idle and enqueue to event queue:
function stopBox(handle)
    global rtbox_info;
    global rtbox_global;
    
    % Box scanning?
    if rtbox_info(handle).boxScanning == 0
        % Box already stopped. Nothing to do:
        return;
    end
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;
    
    % Wait until we can be sure that the box is ready to receive new
    % commands. The deadline is computed so that in the worst conceivable
    % case the box will be ready to receive at least 1 new command byte:
    WaitSecs('UntilTime', rtbox_info(handle).busyUntil);
    
    % Box ready to receive our "disable all event reporting" command code
    % 'a'. Send it blocking:
    [nw, tpost] = IOPort('Write', s, 'a');
    
    % Command submission completed at time 'tpost'. Set a new busyUntil
    % time. After that time, the box should have stopped at the latest:
    rtbox_info(handle).busyUntil = tpost + rtbox_global.maxbusy;
    WaitSecs('UntilTime', rtbox_info(handle).busyUntil);
    
    % Box should be idle now. Submit a marker token which the box can
    % acknowledge, so we know when all receive queues have fully drained.
    % We send all disable codes, as they are no-ops (all events are already
    % disabled) but get acknowledged by the box:
    markerStr = 'dufopa';
    for i=1:length(markerStr)
        IOPort('Write', s, markerStr(i));
        % We delay for 0.1 msecs after each write just to make absolutely
        % sure it gets through:
        WaitSecs(0.0001);
    end
    
    % Marker submitted. After some roundtrip delay, the box will feed this
    % back as the last element in the receive buffer, after all events that
    % were detected and dispatched by the box before the 'a' disable all
    % call took effect.
    
    % Drain the serial receive buffer until the marker token has been
    % received and all queues are completely empty: We allow thist to take
    % as long as neccessary, i.e., absolute timeout is inf-inity, but we
    % abort if more than 10 seconds elapse without any new data. Rationale:
    % If tons of events are pending in the receive buffer, it might take
    % long to fetch all of them, therefore the generous unlimited inf, but
    % if nothing new arrives within 10 seconds then the buffer is clearly
    % empty, no new data is arriving from box and the marker token hasn't
    % been detected (-1), so something went seriously wrong, e.g.,
    % communication failure with box or corrupt marker token received, and
    % we need to recover by aborting the parse-Op:
    parseQueue(handle, -1, inf, 10); 
    
    % Box event scanning and reporting is stopped:
    rtbox_info(handle).boxScanning = 0;
    
    % Ok, the box is completely idle and all serial queues are empty. The
    % software queue contains all remaining events.
    return;
end

% Parse serial port receive queue, dequeue all stored events and store it
% in internal software queue.
function [nadded, tlastadd] = parseQueue(handle, minEvents, timeOut, interEventDelay)
    global rtbox_info;
    
    % Count of added items in this invocation:
    nadded = 0;
    tlastadd = 0;
    
    % Child protection:
    if ~rtbox_info(handle).boxScanning
        warning('PsychRTBox: parseQueue() called but box not in active scanning mode! Driverbug?!?'); %#ok<WNTAG>
    end
    
    % Get handle to serial port:
    s = rtbox_info(handle).handle;

    if minEvents == -1
        untilMarker = 1;
        minEvents = inf;
    else
        untilMarker = 0;
    end
    
    if minEvents == -2
        untilAllAcknowledged = 1;
        minEvents = inf;
    else
        untilAllAcknowledged = 0;
    end
    
    % Setup deadline for stop of parsing, based on timeOut:
    tcurrent = GetSecs;
    timeOut  = tcurrent + timeOut;
    tInterTimeout = tcurrent + interEventDelay;

    % Parse until timeOut reached, or minimum number of requested events
    % dequeued, whatever comes first:
    while (tcurrent <= timeOut) && (minEvents > 0) && (tcurrent <= tInterTimeout)

        % Shall we repeat until all outstanding tokens acknowledged?
        if untilAllAcknowledged && isempty(rtbox_info(handle).ackTokens)
            % Yes. And we're done:
            break;
        end

        % Fetch at most one byte, non-blocking: Also update current
        % timestamp:
        [evid, tcurrent] = IOPort('Read', s, 0, 1);
        
        if isempty(evid)
            % Nothing received. Sleep a msec, slacky, to release cpu:
            tcurrent = WaitSecs('YieldSecs', 0.001);
            continue;
        end
        
        % Something received. Update interevent timeout:
        tInterTimeout = tcurrent + interEventDelay;

        % rtbox_global.eventcodes=[49:2:55 50:2:56 97 48 57 89]; % code for 12 events
        
        % Dispatch depending on type of received data:
        
        % 7-Byte event packet?
        if (evid >= 48 && evid <= 57) || (evid == 89) || (evid == 97)  
            % Read remaining 6 bytes of 7 byte packet, either a timestamp
            % or parameters. We do this blocking as we know there must be 6
            % more bytes pending:
            [b6, tcurrent] = IOPort('Read', s, 1, 6);
            if length(b6)<6
                warning('PsychRTBox: parseQeue: Corrupt 7-byte data packet received from box!'); %#ok<WNTAG>
            else
                % evid 'a' is a special case: An 'a' can be either an
                % acknowledge for disable of all event reporting, or a electronic
                % pulse trigger event :-(
                % Treat this special case of character a == 97, but only if
                % a marker token is expected. Then we might receive the
                % special 'dufopa' signature in b6:
                if untilMarker && (evid == 97) && strcmp(char(b6), 'dufopa')
                    % Marker token detected! This is the end of our
                    % parsing operation, as we know no further data can
                    % be received from the box after this marker.
                    
                    % Btw. the 'dufopa' also corresponds to a specific 6
                    % byte timestamp, so this is ambiguous in theory.
                    % However the timestamp corresponding to 'dufopa' would
                    % correspond to a box uptime (and therefore computer
                    % uptime) of over 30 years -- Unlikely / Impossible to
                    % ever happen in reality, so this is a non-issue.
                    
                    % Break out of parser loop:
                    break;
                end
                
                % If we reach this point, then it is a standard event
                % packet with a b6 6 byte box timestamp.
                
                % Decode timestamp into uncorrected box time in seconds:
                secs=bytes2secs(b6, handle);
                                
                % Buttonstate update needed?
                if (evid >= 49 && evid <= 56)
                    % One of the four pushbuttons: Pressed or released?
                    buttonid = floor((evid - 47)/2);
                    isdown   = mod(evid, 2);

                    % Update live state:
                    rtbox_info(handle).buttons(1, buttonid) = isdown;
                    
                    % Debouncer implementation: Is this button updated
                    % after the debounce deadline which defines the "dead"
                    % interval of ignoring this button?
                    if secs < rtbox_info(handle).buttons(2, buttonid)
                        % This event was received within the "dead
                        % interval" in which state updates for the button
                        % should be ignored for the purpose of button
                        % debouncing. We abort processing for this event,
                        % ie., the event gets discarded:
                        continue;
                    else
                        % Button event received after debounce deadline. We
                        % accept this event and update the debounce
                        % deadline accordingly:
                        rtbox_info(handle).buttons(2, buttonid) = secs + rtbox_info(handle).buttons(3, buttonid);
                    end
                    
                    % Special handling for button press/release events.
                    % These are always reported by box if one of them is
                    % enabled. Check if we should discard some of them:
                    if (isdown && ~rtbox_info(handle).enabled(1)) || (~isdown && ~rtbox_info(handle).enabled(2))
                        % This event shall not be reported to higher
                        % layers. Skip it:
                        continue;
                    end
                end
                
                % Are there any tokens pending for acknowledge?
                if ~isempty(rtbox_info(handle).ackTokens)
                    % Yes. Unacknowledged tokens pending. We discard this
                    % event. We will discard all events until there are no
                    % more unacknowedged tokens around:
                    continue;
                end

                % Enqueue event in internal queue:
                rtbox_info(handle).recQueue(end+1, :) = [ double(evid), secs ];
                
                % Increment total count of added events:
                nadded = nadded + 1;
                
                % Update timestamp of last event add operation:
                tlastadd = GetSecs;
                
                % Decrement count of events to parse:
                minEvents = minEvents - 1;
            end
            
            % Next iteration:
            continue;
        end

        % Any acknowledge tokens to wait for?
        if ~isempty(rtbox_info(handle).ackTokens)
            % Is current event the next expected token?
            if evid == rtbox_info(handle).ackTokens(1)
                % Yes. Dequeue this token.
                rtbox_info(handle).ackTokens = rtbox_info(handle).ackTokens(2:end);
                
                % Next iteration:
                continue;
            end
        end

        % Button state live query result packet? This has the evid 63 == '?'
        % Ok, this doesn't work due to a design-flaw in the protocol. The
        % box sends the data first in the evid byte, then the '?'
        % identifier in the 2nd byte, so our parser can't dispatch
        % reliably. :-(
        %
        % The only solution / workaround is to never submit a '?' query
        % while the box is in scanning mode, so this event can't ever
        % happen. -- We use lot's of ugly magic in buttonQuery() to work
        % around this synchronously...
        %
        % Note: The bug is fixed in firmware v1.4 and later, but that
        % doesn't help us anymore, so we keep the current strategy...
        %
        %
        %         if evid == 63
        %             b1=IOPort('read', s, 1, 1); % ? returns 1 byte with encoded button state.
        %             if isempty(b1) || mod(b1,16)~=0
        %                 warning('PsychRTBox: parseQeue: Corrupt 2-byte data packet received from box!'); %#ok<WNTAG>
        %             else
        %                 % Uppermost 4 bits are button states:
        %                 rtbox_info(handle).buttons(1, :) = [ 0, 0, 0, 0 ];
        %                 % Set corresponding entries for pressed buttons:
        %                 rtbox_info(handle).buttons(1, logical(bitget(b1, 5:8))) = 1;
        %             end
        %
        %             % Decrement count of events to parse:
        %             % This is a non-event in the sense of ther parse function... minEvents = minEvents - 1;
        %
        %             % Next iteration:
        %             continue;
        %         end

        % Next parse iteration:
    end

    if untilAllAcknowledged && ~isempty(rtbox_info(handle).ackTokens)
        warning(['PsychRTBox: parseQeue: Did not receive expected acknowledgements from box! Missing acks: ' char(rtbox_info(handle).ackTokens)]); %#ok<WNTAG>
    end

    % Done.
    return;
end

% Retrieve events from internal driver queue, optionally fetch new items
% from serial port device receive queue if internal queue can't satisfy the
% amount of requested events or usercode wants full coverage of all what is
% there:
%
% - Most of the time you will want to set minItems == MaxItems == nItems
%   with nItems being the exact amount of events you'd like to have within
%   some timeout period (for recovery from sleeping subjects), or nItems ==
%   inf and some reasonable timeouts to get all responses within some
%   defined response interval.
%
% - If the box is not in scanning mode, this will return all data up to
%   'maxItems' events from the internal receive queue instantaneously and
%   the timeouts and minItems don't apply.
%
function [evts, boxtimes] = getEvents(handle, minItems, maxItems, timeOut, maxInterEvent)
    global rtbox_info;

    % Must not operate much longer than tdeadline for given timeOut:
    tcurrent = GetSecs;
    tdeadline = tcurrent + timeOut;
    tInterdeadline = tcurrent + maxInterEvent;
    evAvail   = 0;

    % Box in scanning mode?
    % If box is not in active scanning mode then all serial receive buffers
    % are empty and the box won't produce new data. All we can ever get is
    % what is already in our internal recQueue. In that case, we skip the
    % hardware buffer scan loop.
    if rtbox_info(handle).boxScanning
        % Yes: New data may arrive from box anytime...
        % Repeat fetch & assign operation until timeout, or enough items found:
        while (tcurrent < tdeadline) && (tcurrent < tInterdeadline) && (evAvail < minItems)
            % How many events are stored in our own filtered, debounced event Queue?
            evAvail = size(rtbox_info(handle).recQueue, 1);

            % How many do we need to fetch from box serial receive buffer?
            needMore = minItems - evAvail;

            % Need any?
            if needMore > 0
                % Need to fetch from serial queue:
                if tdeadline - tcurrent < tInterdeadline - tcurrent
                    tscandeadline = tdeadline - tcurrent;
                else
                    tscandeadline = tInterdeadline - tcurrent;
                end
                [nadded, tlastadd] = parseQueue(handle, needMore, tscandeadline, maxInterEvent);
            else
                nadded = 0;
            end

            % Update time:
            tcurrent = GetSecs;

            % Anything added?
            if nadded > 0
                % Yes. Update interevent deadline:
                tInterdeadline = tlastadd + maxInterEvent;
            end

            if (needMore > 0) && (nadded == 0)
                % Nothing received. Sleep a msec, slacky, to release cpu:
                tcurrent = WaitSecs('YieldSecs', 0.001);
            end
            
            % Check again.
        end
    end

    % How many events are stored in our own filtered, debounced event Queue?
    evAvail = size(rtbox_info(handle).recQueue, 1);

    % Hopefully got something. Extract subset of at most maxItems oldest
    % elements:
    if evAvail > maxItems
        evAvail = maxItems;
    end
    
    % Assign first evAvail elements of queue to return arguments:
    if evAvail > 0
        evts     = rtbox_info(handle).recQueue(1:evAvail, 1);
        boxtimes = rtbox_info(handle).recQueue(1:evAvail, 2);

        % Strip them from queue:
        rtbox_info(handle).recQueue = rtbox_info(handle).recQueue(evAvail+1:end, :);
    else
        evts     = [];
        boxtimes = [];
    end
    
    return;
end


% Purge all buffers from stale data:
function purgeRTbox(handle)
    global rtbox_info;

    % Get IOPort handle:
    s = rtbox_info(handle).handle;

    % Stop all processing on box, drain its buffers:
    stopBox(handle);
    
    if IOPort('BytesAvailable', s)
        warning('PsychRTBox: purgeRTBox: Unexpected data in receive buffer after stopBox()! n = %i\n', IOPort('BytesAvailable', s)); %#ok<WNTAG>
    end
    
    % Make additional low-level call to clear all serial port buffers -
    % This is redundant:
    IOPort('Purge', s);
    
    if IOPort('BytesAvailable', s)
        warning('PsychRTBox: purgeRTBox: Unexpected data in receive buffer after serial purge! n = %i\n', IOPort('BytesAvailable', s)); %#ok<WNTAG>
    end

    % Clear the receive event queue of the box:
    rtbox_info(handle).recQueue = [];
    
    % All pending events cleared:
    return;    
end

% Convert 6-byte raw box timestamp x into seconds: boxClockTickIntervalSecs
% is the time unit of the device clock, i.e., the device clock increments
% with a granularity of boxClockTickIntervalSecs:
function secs=bytes2secs(b6, handle)
    global rtbox_info;
    secs = (256.^(5:-1:0) * b6(:) * rtbox_info(handle).boxClockTickIntervalSecs);
end

% Open first available (not yet opened) RT-Box if no specific 'deviceID' is
% given, otherwise try to open box with that 'deviceID'. Store info about
% opened box in slot 'handle' of the rtbox_info struct array:
function openRTBox(deviceID, handle)
    global rtbox_info;
    global rtbox_global; %#ok<NUSED>
    
    % Setup possible port probe list for different operating systems:
    if IsWin
        if ~isempty(strfind(deviceID, 'COM'))
            % COM port name provided: Use specified COM port:
            ports{1} = deviceID;
        else
            % Enumerate all possible ports as candidates for box:
            % Suppose RTBox not assigned to COM1 or COM2, as these are
            % usually the native serial ports, if any:
            ports=cellstr(num2str((3:256)','\\\\.\\COM%i'));
            ports=strtrim(ports); % needed for Matlab R2009b and later.
        end
        
        nPorts=length(ports);
    end

    if IsOSX || IsLinux
        if ~isempty(strfind(deviceID, '/dev/'))
            % Explicit port name provided: Use it.
            ports{1} = deviceID;
            nPorts = 1;
        else
            % Enumerate all available USB-Serial ports as candidates for box:
            if IsOSX
                candidateports=dir('/dev/cu.usbserial*');
            else
                candidateports=dir('/dev/ttyUSB*');
            end

            for i=1:length(candidateports)
                ports{i} = ['/dev/', candidateports(i).name];
            end

            nPorts = length(candidateports);
        end
        
    end

    % Any found?
    if nPorts==0
        RTboxError('noUSBserial');
    end

    % Silence IOPorts output to prevent visual clutter:
    verbo=IOPort('Verbosity', 0);

    % Prepare device detection loop:
    deviceFound=0;
    rec=struct('avail','','busy',''); % for error record only
    trycount = 0;
    
    % Retry device detection procedure up to a total of 2 times:
    % Sometimes a freshly connected RTBox has a little hickup at first try,
    % but it reliably works at the 2nd iteration:
    while (trycount < 2) && (~deviceFound)    
        % Scan each possible port if it is the response box:
        for ic=1:nPorts
            % Device file name of port to probe:
            port = ports{ic};

            % Already opened?
            if strmatch(port, {rtbox_info(1:handle-1).portname},'exact')
                % Yes. Skip this candidate:
                continue;
            end
            
            if IsLinux
                % Switch serial port to low-latency mode on Linux:
                % This is extra-paranoid, as our udev rules and the
                % 'ReceiveLatency=0.0001' setting below does already the same.
                system(sprintf('setserial %s low_latency', port));
                WaitSecs('YieldSecs', 0.100);
            end

            % Try to open port: We open at maximum supported baud rate of 115200
            % bits, use a timeout for blocking reads of 1 second, and request
            % low-latency polling for write completion if the IOPort('Write')
            % command uses a polling method for waiting for write completion. Set
            % the "sleep time" between consecutive polls to 0.0001 seconds = 0.1
            % msecs. That is good enough for our purpose and still prevents system
            % overload on OS/X and Linux. We also set the ReceiveLatency to 0.1 msecs,
            % aka 0.0001 seconds: This parameter is silently ignored on Windows, honored
            % in some way by some serial port drivers on OS/X (well, maybe, who knows?),
            % and on Linux any value <= 1 msecs enables ASYNC_LOW_LATENCY mode on serial
            % ports, ie. some low-latency optimizations. E.g., with the FTDI chips used in
            % RTBox it will automatically set the chips latency timer to its minimum of 1 msec:
            [s errmsg]=IOPort('OpenSerialPort', port, 'BaudRate=115200 ReceiveTimeout=1.0 PollLatency=0.0001 ReceiveLatency=0.0001 ');

            % Worked?
            if s>=0
                % Device open succeeded. Test if it is really our box and not some
                % other serial-USB device:

                % Flush all transmit/receive queues:
                % N.B. This flushes receive and send buffers on Windows and OS/X,
                % but doesn't flush the send buffers on Linux -- not implemented at
                % the FTDI driver level as of Linux 2.6.22:
                IOPort('Purge', s);

                % Read out whatever junk maybe in the input buffer:
                IOPort('Read', s, 0);

                % Write the 'X' command code to ask box for its identity and to
                % switch it into advanced mode with the clever protocol we
                % use. The box powers up in "simple" mode, where most cool
                % functionality and protocol is disabled but it acts in a
                % way that is compatible to E-Prime et al. -- Just sends 1
                % Byte event code whenever some button or input port
                % changes state:
                IOPort('Write', s, 'X');

                % Wait blocking with 1 sec timeout (see above) for id string response from box:
                idn=char(IOPort('Read', s, 1, 21));
                
                % Expected response is a 21-Bytes string of format 'USTCRTBOX,115200,v?.?'
                % with ?.? being the major and minor firmware/box revision.
                if strfind(idn,'USTCRTBOX')
                    % Found device:
                    deviceFound=1;
                    break;
                end

                % Not our box :-(
                IOPort('Close', s);

                % Store port as existent but not used by us:
                rec.avail{end+1}=port; %#ok

            elseif isempty(strfind(errmsg,'ENOENT'))
                % Failed to open port, but not with error code ENOENT.
                
                % Open failed, but port exists. That means it is busy - used by
                % ourselves or some other process:
                rec.busy{end+1}=port; %#ok

                if isempty(strfind(errmsg,'EPERM')) && isempty(strfind(errmsg,'EBUSY'))
                    % No permission or busy error. Something else went
                    % wrong. Tell user:
                    RTboxWarn('openError', errmsg);
                end
            end
            % Scan next candidate port:
        end
       
        % Tried trycount times to perform open operation:
        trycount = trycount + 1;
        
        % Retry, if neccessary:
    end
    
    % Done with scan. Restore normal level of debug output of IOPort, as
    % selected by external usercode or defaults:
    IOPort('Verbosity',verbo);

    % Found a suitable RTBox within number of retries?
    if ~deviceFound
        % Nope. Bail out:
        RTboxError('noDevice', rec, rtbox_info);
    end

    % Ok, found a box and opened connection. Setup its device info structure:
    
    % First the default settings...
    rtbox_info(handle)=struct('events',{{'1' '2' '3' '4' '1up' '2up' '3up' '4up' 'pulse' 'light' 'tr' 'serial'}},...
                              'enabled',[], 'ID','','handle',-1,'portname',[],'sync',[],'version',[],'clkRatio',1,'verbosity',3, ...
                              'busyUntil', 0, 'boxScanning', 0, 'ackTokens', [], 'buttons', [0 0 0 0; 0 0 0 0; 0.05 0.05 0.05 0.05], ...
                              'syncSamples', [], 'recQueue', [], 'boxClockTickIntervalSecs', 1/115200);

    % Enabled events at start:
    rtbox_info(handle).enabled=logical([1 0 0 0 0 0]);
    % IOPort serial port handle:
    rtbox_info(handle).handle=s;
    % DeviceID:
    rtbox_info(handle).ID=deviceID; % store handle and deviceID
    % Portname:
    rtbox_info(handle).portname=port;
    % Box version:
    ind = strfind(idn, ',v') + 2;
    rtbox_info(handle).version = str2double(char(idn(ind + (0:2))));
    
    % Duration of a single clock tick of the RTBox microprocessor clock:
    % Defaults to 1/115200 second, but will be overridden to 1/921600 second
    % if the spec string is found as part of the idn string. This finer
    % resolution is present in the latest RTBox devices (Firmware v1.4:
    if strfind(idn,'921600')
        rtbox_info(handle).boxClockTickIntervalSecs = 1/921600;
    end
    
    % Init clock-ratio to an uncalibrated 1.0:
    rtbox_info(handle).clkRatio=1;
    
    % Start with all scanning disabled. We do the initial disable via
    % low-level writes to make absolutely sure that the box is "off", as
    % all other utility routines rely on a well defined initial off state.

    % Read out whatever junk maybe in the input buffer:
    IOPort('Read', s, 0);
    
    % Submit disable all command:
    if IOPort('Write', s, 'a')~=1
        RTBoxError('[In event disable]: Failed to initialize RTBox to proper startup settings!');
    end
    
    % Wait 200 msecs:
    WaitSecs(0.2);
    
    % Go into an infinite loop that repeats until nothing available in
    % receive buffer:
    while 1
        curItem = char(IOPort('Read', s, 1, 1));
        if IOPort('BytesAvailable', s) == 0
            break;
        end
        
        WaitSecs(0.050);
    end

    % Last item == disable all acknowledge?
    if curItem ~= 'a'
        RTBoxError('[In event disable II]: Failed to initialize RTBox to proper startup settings!');
    end
    
    % Ok, now the box is fully idle and we return with all event scanning
    % disabled.
    
    % Device open and initialized.
    fprintf('PsychRTBox: RTBox device "%s" opened on serial port device %s.\n', deviceID, port);
end

% put verbose error message here, to make main code cleaner
function RTboxError(err,varargin)
    global rtbox_info;
    global rtbox_global;

    switch err
        case 'invalidTTL'
            str='TTL output must be integer from 0 to 15, or bits from ''0000'' to ''1111''.';
        case 'noUSBserial'
            str=' No connected USB-Serial devices found. Is your device connected and the USB-Serial driver properly installed?';
        case 'noDevice'
            p=varargin{1};
            str=sprintf(' No RTBox on available ports found or no USB-Serial port found at all, or no port with the given name exists.');

            if ~isempty(p.avail) % have available ports
                str=sprintf([' USB-Serial Port(s) available: %s, but failed to identify any of them as an RTBox device.\n' ...
                    ' Is any of them the RTBox device? If yes, try again.\n' ...
                    ' It may help to unplug, then replug the device.\n\n'],cell2str(p.avail));
            end

            if ~isempty(p.busy) % have busy ports
                str=sprintf(['%s USB-Serial Port(s) unavailable:\n%s\nProbably already in use by another program.\n' ...
                    'Is any of the given busy ports the RTBox device? If yes, try ''clear all'' to close the ports and then retry.\n'], str, cell2str(p.busy));
            end

            % Have opened RTbox?
            if rtbox_global.nrOpen > 0
                str=sprintf('%s\nAlready opened RTBox device(s): ', str);
                for i=1:length(rtbox_info)
                    str=sprintf('%s %s at %s,\n',str, char(rtbox_info(i).ID), char(rtbox_info(i).portname));
                end
                str(end)='.';
            end
        case 'unknownCmd'
            str=sprintf(['Unknown subcommand or trigger: ''%s''.\n'...
                'The first string input must be one of the commands or events:\n%s,\n%s.'],...
                varargin{1},cell2str(varargin{2}),cell2str(varargin{3}));
        case 'invalidButtonNames'
            str=sprintf('ButtonNames requires a cellstr containing four button names.');
        case 'invalidButtonName'
            str=sprintf('Invalid button name: %s.',varargin{1});
        case 'notRespond'
            str=sprintf('Failed to communicate with the device due to communication error. Try to close and re-connect the device.');
        case 'invalidEnable'
            str=sprintf('Valid events for enable/disable: %s.',cell2str(varargin{1}));
        case 'triggerDisabled'
            str=sprintf('Trigger is not enabled. You need to enable ''%s''.',varargin{1});
        otherwise, str=err;
    end
    
    error(sprintf('PsychRTBox:%s %s', err, WrapString(str) )); %#ok<SPERR>
end

% Show warning message, but code will keep running
function RTboxWarn(err,varargin)
    switch err
        case 'invalidTTLwidth'
            str=sprintf('Supported TTL width is from %.2fe-3 to %.2fe-3 s .',[1 255]*varargin{1}*1000);
        case 'widthOffset'
            str=sprintf('TTL width will be about %.5f s', varargin{1});
        case 'invalidScanNum'
            str=sprintf('The scan number for HardwareDebounce must be an integer from 1 to 255.');
        case 'notSupported'
            str=sprintf('The command %s is supported only for Firmware version V%.1f or later.',varargin{1:2});
        case 'invalidEvent'
            str=sprintf('%4i',varargin{1});
            str=sprintf('PsychRTBox: Events not recognized: %s.\nPlease do PsychRTBox(''clear'') before showing stimulus.\nGetSecs = %.1f',str,GetSecs);
        case 'noTrigger'
            str=sprintf('PsychRTBox: Trigger ''%s'' not detected. GetSecs = %.1f', varargin{1}, GetSecs);
        case 'poorSync'
            str=sprintf('PsychRTBox: Unreliable synchronization detected (variation range %.1f ms).\nGetSecs = %.1f', varargin{1}*1000,GetSecs);
        case 'allDisabled'
            str=sprintf('PsychRTBox: All event detection has been disabled for %s.', varargin{1});
        case 'openError'
            str=sprintf('PsychRTBox: Some unknown error occured while trying to open port. The IOPOrt driver says: %s', varargin{1});            
        otherwise
            str=sprintf('PsychRTBox: %s. GetSecs = %.1f',err,GetSecs);
    end

    str=WrapString(str);
    warning(['PsychRTBox:' err], str);
end

% return str from cellstr for printing, also remove port path
function str=cell2str(Cstr)
    if isempty(Cstr)
        str='';
        return;
    end

    str='';
    for i=1:length(Cstr)
        str = [str Cstr{i} ' ']; %#ok<AGROW>
    end
end
