function varargout = Eyelink(varargin)
%  The EyelinkToolbox can be used to ceate eye-movement experiments and
%  control the SR-Research Eyelink gazetrackers (https://www.sr-research.com/)
%  from within Octave and Matlab.
%
%  It is incorporated into the PsychToolbox (http://www.psychtoolbox.org/)
%  and uses the functions provided by the PsychToolbox for graphics and sound.
%
%  To use Eyelink, you must also install the SR-Research runtime libraries,
%  downloadable from the Support Section of their website (see above).
%
%  For a complete list of available functions type "Eyelink" in the command
%  window. For an explanation of any particular Eyelink function just add a
%  question mark "?" after a command.
%
%  E.g. for 'Initialize', try either of these equivalent forms:
%  Eyelink('Initialize?')
%  Eyelink Initialize?
%
%  [optional arguments]:
%  Brackets in the function list, e.g. [remport], indicate optional arguments,
%  not matrices. Optional arguments must be in order, without omitting earlier
%  ones.
%
%  If you need examples to get you started, check out the EyelinkDemos
%  folder (help EyelinkDemos).
%
%  More information on this toolbox can be found in the file:
%  EyelinkToolbox/contents.m  (help EyelinkToolbox)
%
%  The EyeLink Toolbox requires the EyeLink Developers Kit installed to
%  provide the dependency libraries. Please register a user account (free)
%  at https://www.sr-support.com for access to all support related
%  information and software downlaods from SR-Research
%

AssertMex('Eyelink.m');
