%  EyelinkToolbox.
%
%  The EyelinkToolbox can be used to ceate eye-movement experiments and
%  control the SR-Research Eyelink(c) gazetrackers
%  (http://www.eyelinkinfo.com/) from within Matlab and Octave.
%
%  It is incorporated into the PsychToolbox (http://www.psychtoolbox.org/).
%
%  It provides an interface to the Eyelink Gazetracker.
%
%  It uses a very similar  approach as the PsychToolBox and uses the
%  functions provided by the PsychToolBox for graphics and sound.
%
%  How to start using it:
%
%  The EyelinkToolbox is a collection of m-files and a Mex file.
%  The EyelinkToolbox core is the mex/dll-file called Eyelink. The
%  EyelinkToolbox uses the same approach as the SCREEN mex/dll function
%  provided in the PsychToolBox. Most of the Eyelink functions that are
%  avaible in C are also available under Matlab. In addition, the
%  EyelinkToolbox provides a number of wrapper functions to simplify
%  creating an eye-tracking program. The main functionality of the 
%  EyelinkToolbox is demonstrated in a number of demos. This is probably
%  where you want to start if you are new to the EyelinkToolbox.
%
%  To use Eyelink, you must also install the SR-Research runtime libraries,
%  downloadable from the Support Section of their website (see above).
%
%  For a complete list of available functions type "Eyelink" in the
%  command window. For an explanation of any particular eyelink function
%  just add a question mark "?" after a command. E.g. for help on the 
%  function 'Initialize', try either of these equivalent forms:
%
%  Eyelink('Initialize?')
%  Eyelink Initialize?
%
%  Contents of the EyelinkToolbox folder:
%
%  -Contents.m: this file.
%
%  -Changes.m: file that documents changes to the EyelinkToolbox
%
%  -EyelinkBasic
%  Eyelink mex files and a collection of m-files that provide the
%  core functionality of the toolbox.
%
%  -EyelinkDemos
%  EyelinkShortDemos is the best place to start. Contains relatively 
%  simple demos. The SR-Research Demo folder is a slightly more
%  elaborate example. The folder BRMIC contains the original example
%  as shown in the EyelinkToolbox article in BRMIC (Cornelissen,
%  Peters and Palmer 2002). It is obsolete in many ways but left in
%  the toolbox for historic reasons.
%
%  -EyelinkOneLiners
%  Handy routines that provide less essential functionality
%
%  We welcome any observations, suggestions that may help us improve this
%  toolbox. If any particular function you need is still missing, let us know
%  and we'll try to incorporate it into a next release. In case you decide
%  to add or modify functions yourself, please send us the modified code.
%  If you think you've found a bug, please tell us. It will help greatly
%  if you can supply a  minimal-length program that exhibits the bug.
%  If you're interested in the c-source code of this project, it is part
%  of the source code of the PsychToolbox (http://www.psychtoolbox.org/).
%
%  ******Citing the use of the Eyelink Toolbox******
%
%  We would very much appreciate it if you would (also) cite the EyelinkToolbox
%  paper when you have used the toolbox in your work. You could write something
%  like:
%
%  We wrote our experiments in Matlab, using the Psychophysics and Eyelink
%  Toolbox extensions (Kleiner 2007; Brainard, 1997; Pelli, 1997; Cornelissen,
%  Peters & Palmer, 2002; see http://psychtoolbox.org).
%
%  Kleiner M, Brainard D, Pelli D, 2007, "What's new in Psychtoolbox-3?"
%  Perception 36 ECVP Abstract Supplement.
%
%  Brainard, D. H. (1997) The Psychophysics Toolbox, Spatial Vision 10:433-436.
%  Cornelissen, F.W., Peters. E., Palmer, J. (2002).
%
%  The Eyelink Toolbox: Eye tracking with MATLAB and the Psychophysics
%  Toolbox. Behavior Research Methods, Instruments & Computers, 34,
%  613-617.
%
%  Pelli, D. G. (1997) The VideoToolbox software for visual psychophysics:
%  Transforming numbers into movies, Spatial Vision 10:437-442.
%  This toolbox is also actively supported and co-developed by SR-Research
%  (http://www.eyelinkinfo.com), the manufacturer of the Eyelink
%  gazetrackers. You can contact them about it. However, if your problem
%  is likely of interest to other users as well, then please post
%  questions on the PsychToolbox mailinglist, which is also monitored by
%  SR-Research.
%
%  The EyelinkToolbox and its source code are copyright 
%  Enno Peters, Frans Cornelissen and John Palmer. Additional copyrights
%  for some contributed functions are SR-Research and the respective
%  individual developers of those functions.
%
%  (Type "Eyelink" or Eyelink('Version') for authors). The
%  toolbox is freely useable and redistributable under the terms of the
%  MIT license. The exact license text is included inside License.txt in
%  the root folder of any Psychtoolbox-3 installation.
%
%  Disclaimer: we cannot be hold responsible for any damage that may
%  (appear to) be caused by the use of this toolbox. Use at
%  your own risk!
%
%
%  Frans W. Cornelissen
%  22nd July 2010
%  email: f.w.cornelissen@med.umcg.nl
%
%  Enno Peters, Frans Cornelissen and John Palmer
%  Groningen, 27-11-2002
%
