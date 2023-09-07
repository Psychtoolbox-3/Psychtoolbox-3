% EyelinkToolbox:EyelinkDemos:SR-ResearchDemos:EyeLink_MRI_BlockRecord
% 6 trials are presented in 2 blocks of 3 trials. Trial duration is 5.5s during which a 4s stimulus is presented.
% A block starts with a drift-check followed by presentation of central crosshairs. Eye movements are recorded while 
% waiting for an MRI trigger (keyboard key 't' in this demo). The stimulus is presented when trigger is received.
% A fixed ITI is maintained by presenting crosshairs between each 4s stimulus. Eye movements are recorded throughout
% an entire block rather than on a trial-by-trial basis. 
%
% In STEP 5 it is shown how to:
% - shrink the spread of the calibration/validation targets so they are all visible if the MRI bore blocks part of the screen
% - apply an optional online drift correction (see EyeLink 1000 Plus User Manual section 3.11.2)
%
% Usage:
% Eyelink_MRI_BlockRecord(screenNumber)