% EyelinkToolbox:EyelinkDemos:SR-ResearchDemos
%
% EyeLink systems allow for rich integration between a Display PC (the computer running Psychtoolbox) and the EyeLink Host PC via an ethernet connection using
% the EyeLink application programming interface (API). The API can be installed by downloading the EyeLink Developers Kit from our Support site:
%
% <a href="https://www.sr-research.com/support/thread-13.html">EyeLink Developers Kit / API Downloads (Windows, macOS, Linux)</a>
%
% You will also need to configure your Display PC network settings. See this link for more details: 
%
% <a href="https://www.sr-research.com/support/thread-58.html">How to configure a Display PC</a>
%
% Psychtoolbox interfaces with the EyeLink API via the 'Eyelink Toolbox' which can be found in the Psychtoolbox root folder > PsychHardware.
% You can add EyeLink integration code to your own Psychtoolbox script to interact with the Host PC so your experiment can:
%    - initialize an EyeLink connection
%    - open an EyeLink data file (edf) on the Host PC and name it
%    - set some EyeLink-specific parameters /change various default options
%    - put Host PC in 'Camera Setup' mode and allow for:
%        - transfer of eye image to Display PC for ease of participant setup 
%        - presentation of targets on Display PC monitor for participant calibration / validation
%    - transfer text, image and/or other graphics to the Host PC at the beginning of each trial for experimenter feedback
%    - do a drift-check / correction
%    - start recording eye movement data at the beginning of each trial (or block)
%    - during a trial write messages in the edf file to mark time of events
%    - write messages in the edf file for rich integration with Data Viewer
%    - stop recording at the end of each trial (or block)
%    - towards the end of a session close the edf and transfer a copy to the Display PC
%    - close the EyeLink connection
%
% EyeLink integration also allows for online access of eye movement data for gaze-contingent experiments.
%
% In a typical task eye movements are recorded onto the EyeLink Host PC on a trial-by-trial basis: recording starts at the beginning of a trial 
% and ends at the end of a trial. This allows for an optional drift-check/correction and for transferring of images to the Host PC between trials. 
% However this trial-based recording might not be suitable for paradigms that require a fixed ITI. In such cases it is possible to record 
% continuously throughout a block or session.  
%
% The following is a list of Psychtoolbox demos with EyeLink integration that are included with the EyeLink Toolbox:
%
% SimplePicture - Eyelink_SimplePicture.m
%
%     A simple EyeLink integration demo that records eye movements while an image is presented on the screen. Each trial ends when the 
%     space bar or a button is pressed.
%
%     Illustrates how to: 
%        - use a non-default IP address for the Host PC
%        - replace the default calibration/validation/drift correct bull's eye target with an image (target recommended by Thaler et al., 2013)
%        - set calibration target and feedback beep sounds on / off
%        - change edf file name when a copy is transferred to the Display PC
%
% GazeContingent - Accessing sample and event gaze data online in real-time over the Ethernet link
%
% Basic EyeLink integration involves saving eye movement data onto the EyeLink Host PC and transferring a copy of the EyeLink data file to the Display PC
% at the end of the session.
% EyeLink eye movement data consists of SAMPLES (eye position recorded at every sample - based on tracker sampling rate)
% and EVENTS (saccades, fixations, blinks detected by the EyeLink online parser) 
% However it is also possible to have Psychtoolbox access eye movement data online for gaze-contingent tasks while still saving data onto the Host PC
%
% There are two ways of accessing eye data online: 
%
% (1) Accessing buffered data by looping through the function pair Eyelink('GetNextDataType') / Eyelink('GetFloatData') to access buffered EVENTS and SAMPLES
%    Use buffered data if you need to:
%    a) grab every single consecutive SAMPLE online without missing any
%    b) grab EVENT data (fixation/saccade/blink events) online
%
%    Note that buffered event data can take some time to be available online due to the time involved
%    in calculating velocity/acceleration. If you need to retrieve online gaze
%    position as fast as possible and/or you don't need to get all available SAMPLES or other
%    EVENTS, then use option (2)
%
%    When using buffered data it is advisable to initially loop through the function pair for ~100ms to clear
%    old data from the buffer before using it, thus allowing access to the most recent EVENTS or SAMPLES
%
%    Demos using buffered data:
%
%    BufferedEndSacEvents - Eyelink_BufferedEndSacEvents.m    
%        A simple EyeLink gaze-contingent demo showing how to retrieve online EVENTS from a buffer.
%        In each trial an image is presented with a red gaze-contingent dot overlaid on top.
%        The dot's location is updated online based on the end x y coordinates of each saccade detected.
%        Each trial ends when the space bar is pressed
%
%    BufferedFixUpdateEvents - Eyelink_BufferedFixUpdateEvents.m
%        A simple EyeLink gaze-contingent demo showing how to retrieve online events from a buffer.
%        In each trial an image is presented with a red gaze-contingent dot overlaid on top.
%        The dot's location is based on the average x y coordinates of fixations updated online every 50ms via a FIXUPDATE event.
%        See EyeLink Programmers Guide manual > Experiment Templates Overview > Control > Fixation Update Events
%        Each trial ends when the space bar is pressed.
%
%        FIXUPDATE events send updates about a current fixation at regular intervals. By default an interval of 50ms is used. 
%        The first update is sent one update interval after the start of the fixation, and the last is sent at
%        the end of the fixation. This demo uses FIXUPDATE events to get the averaged gaze x y position across each fixation interval.        
%
%    FixWindowBufferedSamples - Eyelink_FixWindowBufferedSamples.m
%        EyeLink gaze-contingent demo that shows how to retrieve online gaze SAMPLES from a buffer.
%        In each trial central crosshairs are shown until gaze is detected continuously within a central 
%        square window for 500ms or until the space bar is pressed. An image is 
%        then presented until the space bar is pressed to end the trial.
%
% (2) Fast access of SAMPLES by looping through the function pair Eyelink('NewFloatSampleAvailable') / Eyelink('NewestFloatSample') to access the most 
%    recent SAMPLES online.
%    Use this option if you need to retrieve online eye position (e.g. GAZE data) as fast as possible.
%    This option may not necessarily retrieve every single consecutive SAMPLE (this depends on how fast your loop is executed by the Display PC) 
%    and this option does not allow for EVENT checking.
%
%    Demos using fast samples (option 2 above):          
%
%    GCfastSamples - Eyelink_GCfastSamples.m
%        A simple EyeLink gaze-contingent demo showing how to retrieve fast online SAMPLES.
%        In each trial an image is presented with a red gaze-contingent dot overlaid on top.
%        The dot's location is updated online based on the x y coordinates of the latest gaze SAMPLE retrieved online.
%        Each trial ends when the space bar is pressed.
%
%    FixWindowFastSamples - Eyelink_FixWindowFastSamples.m
%        EyeLink gaze-contingent demo showing how to retrieve fast gaze samples online.
%        In each trial central crosshairs are shown until gaze is detected continuously within a central 
%        square window for 500ms or until the space bar is pressed. An image is 
%        then presented until the space bar is pressed to end the trial.
%
%
% SimpleVideo - Eyelink_SimpleVideo.m
%
%     Simple video demo with EyeLink integration and animated calibration / drift-check/correction targets.
%     In each trial eye movements are recorded while a video stimulus is presented on the screen.
%     Each trial ends when the space bar is pressed or video stops playing. A different drift-check/correction
%     animated target is used in each of the 3 trials.
%
%     Illustrates how to send messages to allow for a video file to be played back in Data Viewer's 'Trial Play Back Animation' view.
%     Shows how to use animated targets: 
%        - replacing the default calibration/validation/drift-check (or drift-correction) targets with a video file
%        - updating the drift-check/correction video file on a trial-by-trial basis
%
%
% StereoPicture - Eyelink_StereoPicture.m
%
%     EyeLink integration demo for stereo presentation.
%     Records eye movements passively while presenting a stereo stimulus. Supports both split-screen mode 
%     and dual-monitor setup.
%     Each trial ends when the space bar is pressed.
%     Data Viewer integration with both left and right eyes superimposed on the same eye window view.
%     Illustrates how to calibrate on both sides of a split-screen or both monitors in dual monitor setup.
%
% MRI_BlockRecord - Eyelink_MRI_BlockRecord.m
%
%     Simple MRI demo with EyeLink integration.
%     6 trials are presented in 2 blocks of 3 trials. Trial duration is 5.5s during which a 4s stimulus is presented.
%     A block starts with a drift-check followed by presentation of central crosshairs. Eye movements are recorded while 
%     waiting for an MRI trigger (keyboard key 't' in this demo). The stimulus is presented when the trigger is received.
%     A fixed ITI is maintained by presenting crosshairs between each 4s stimulus. Eye movements are recorded throughout
%     an entire block rather than on a trial-by-trial basis. 
%
%     Illustrates how to:
%        - shrink the spread of the calibration/validation targets so they are all visible if the MRI bore blocks part of the screen
%        - apply an optional online drift correction (see EyeLink 1000 Plus User Manual section 3.11.2)
%
% PursuitTarget - Eyelink_PursuitTarget.m
%
%     A smooth pursuit EyeLink integration demo that records eye movements 
%     while a target moves sinusoidally across the screen. Each trial ends after 5s.
%
%     Illustrates how to:
%        - change the drift-check/correction target location before each trial
%        - create a moving target for Data Viewer's Play Back Animation view
%        - create dynamic target location for Data Viewer's Temporal Graph view and sample reports
%        - create target dynamic interest areas for Data Viewer
%
%
% Thaler L, Schutz AC, Goodale MA, Gegenfurtner KR. What is the best fixation target? The effect of target shape on stability of fixational eye movements. Vision Res. 2013; 76: 31–42.
%