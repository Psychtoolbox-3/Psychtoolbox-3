% Psychtoolbox.
% Version 3.0.13      1 June 2016
%
% The Psychophysics Toolbox is a collection of Octave and Matlab functions
% that extend the language to give you exquisite control of your computer
% to test observers with research-grade dynamic stimuli.
%
% Every Psychtoolbox function has its own documentation available through
% the Matlab HELP command, and, in the case of MEX files, through the
% function itself. We've tried hard to make these individual bits of
% documentation clear, accurate, and complete. We're grateful to receive
% corrections.
%
% We've mostly followed Mathworks's help-text conventions, but note that in
% our help text we designate optional arguments to function calls by
% embracing them with square brackets. You're not meant to include these
% brackets when you actually call the function. For example, "help Snd"
% will tell you this:
%
%      err = Snd(command,[sig],[rate])
%
% meaning that the "command" argument is required and the "sig" and "rate"
% arguments are optional. A typical call to Snd looks like this:
%
%      Snd('Play', mysoundsignal)
%
% The most important and generally useful functions, especially Screen, are
% in PsychBasic. The PsychDemos will show you how to do useful things in
% Matlab with the Psychtoolbox. Try these:
%
%      help PsychDemos       % Triple-click me & hit enter.
%      help PsychBasic
%      help Screen
%      Screen OpenWindow?
%
% The folder PsychDocumentation contains bits of documentation about
% technical implementation details, how to solve specific tasks and how to
% troubleshoot common problems.
%
% To read more, look at the "intro" and "docs" pages at the web site. The
% web site also has advice about getting help, the Psychtoolbox forum, and
% the latest information about bugs and new releases:
%
% web http://psychtoolbox.org/; % triple-click me & hit enter
%
% Peter Scarfe has written a large set of nice and beginner friendly
% Psychtoolbox tutorial scripts.
%
% You can find them at http://peterscarfe.com/ptbtutorials.html
%
% If you want to acknowledge use of this software when you publish your
% research, you might say something like this, "We wrote our experiments in
% Matlab, using the Psychophysics Toolbox extensions (Kleiner et al, 2007)"
%
% Kleiner M, Brainard D, Pelli D, 2007, "What's new in Psychtoolbox-3?"
% Perception 36 ECVP Abstract Supplement.
%
% If you find this software useful then please consider donating some money
% to support its ongoing maintenance and development. See:
%
% http://psychtoolbox.org/donations
%
%
% Folders and their contents:
%
% PsychAlpha             - Under development. Experimental, risky, undocumented.
% PsychAlphaBlending     - OpenGL alpha-channel blending utilities and constants.  
% PsychBasic             - Basic support routines for psychophysics.
% PsychCal               - Calibrate your video monitors.
% PsychCalDemoData       - Demo calibration data.
% PsychColorimetric      - Colorimetric calculations.
% PsychColorimetricData  - Standard colorimetric data.
% PsychContributed       - Contributed programs.
% PsychDemos             - Show how to use the Psychtoolbox.
% PsychDocumentation     - Documentation about specific topics.
% PsychFiles             - Process text files.
% PsychGamma             - Fit monitor gamma functions.
% PsychGLImageProcessing - Built-in image processing via graphics hardware.
% PsychGPGPU             - General purpose accelerated computing on GPUs.
% PsychHardware          - Interface to plug-in hardware.
% PsychInitialize        - Initialize and deinitialize MATLAB
% PsychMatlabTests       - Document the few bugs in Matlab 5.2.1.
% PsychObsolete          - Obsolete routines, still present for compatibility.
% PsychOpenGL            - Routines for low-level access to OpenGL 3D graphics.
% Psychometric           - Psychometric function fitting.
% PsychOneliners         - Trivial, but handy, functions.
% PsychOptics            - Optics calculations, mostly for human optics.
% PsychPriority          - Priority and Rush.  (formerly within PsychBasic)
% PsychProbability       - Probability and statistics.
% PsychRadiometric       - Radiometric/photometric calculations and conversions.
% PsychRects             - Manipulate rectangles for drawing.
% PsychSignal            - Signal processing and math routines.
% PsychStairCase         - Adaptive staircase procedure.
% PsychTests             - Evaluate performance of software and hardware.
% PsychVideoCapture      - Functions for realtime video capture.
% Quest                  - Threshold estimation procedure.


% The OS 9 version was this:
%   PsychBasic            - Basic support routines for psychophysics.
%   PsychBeta             - Under development.  May be useful but no promises.
%   PsychCal              - Calibrate your video monitors.
%   PsychCalDemoData      - Demo calibration data.
%   PsychClipboard        - Extend interface between Matlab and clipboard.
%   PsychColorimetric     - Colorimetric calculations.
%   PsychColorimetricData - Standard colorimetric data.
%   PsychContributed      - Contributed programs, not supported by us.
%   PsychDemos            - Show how to use the Psychtoolbox.
%   PsychFiles            - Process text files.
%   PsychGamma            - Fit monitor gamma functions.
%   PsychHardware         - Interface to optional hardware.
%   PsychMatlabTests      - Document the few bugs in Matlab 5.2.1.
%   PsychObsolete         - Obsolete routines, still present for compatibility.
%   Psychometric          - Psychometric function fitting.
%   PsychOneliners        - Trivial, but handy, functions.
%   PsychOptics           - Optics calculations, mostly for human optics.
%   PsychProbability      - Probability and statistics.
%   PsychQuest            - Threshold estimation procedure.
%   PsychRects            - Manipulate rectangles for drawing.
%   PsychSignal           - Signal processing and math routines.
%   PsychTests            - Evaluate performance of software and hardware.


