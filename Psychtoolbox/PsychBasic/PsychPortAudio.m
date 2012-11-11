function varargout = PsychPortAudio(varargin)
% PsychPortAudio - High precision sound driver for Psychtoolbox-3.
%
% PsychPortAudio is a new sounddriver for PTB-3. It is meant to become a
% replacement for all other Matlab based sound drivers and PTB's old SND()
% function.
%
% PsychPortAudio provides the following features:
%
% - Allows instant start of sound playback with a very low onset latency
%   compared to other sound drivers (on well working hardware).
%
% - Allows start of playback at a scheduled future system time: E.g.,
%   schedule sound onset for a specific time in the future (e.g., visual
%   stimulus onset time), then do other things in your Matlab code.
%   Scheduled start of playback can be accurate to the sub-millisecond level
%   on some system setups.
%
% - Wait for sound onset, or continue with execution of your code
%   immediately.
%
% - Asynchronous operation: Sound playback works in the background while
%   your code continues to do other things.
%
% - Infinitely repeating playback, or playback of a sound for 'n' times.
%
% - Returns timestamps and status for all crucial events.
%
% - Support multi-channel devices, e.g., 8-channel sound cards.
%
% - Supports multi-channel sound capture and full-duplex capture
%   and playback of sound on some systems.
%
% - Enumerate, open and use multiple sound cards in parallel.
%
% - Reliable (compared to Matlabs sound facilities).
%
% - Efficient, causes only very low cpu load.
%
% See the "help InitializePsychSound" for more info on low-latency
% configurations. See "help BasicSoundOutputDemo" for a very basic demo of
% sound output (without special emphasis on low-latency). See
% "BasicSoundInputDemo" for a basic demo of sound capture.
% "BasicSoundFeedbackDemo" shows how to implement a simple audio feedback
% loop with controllable delay.
%
% "PsychPortAudioTimingTest" is a script that we used for testing PA's
% sound onset latency and accuracy. It also serves as an example on how to
% get perfectly synched audio-visual stimulus onsets.
%
% Type "PsychPortAudio" for an overview of supported subfunctions and
% "PsychPortAudio Subfunctionname?" for help on a specific subfunction.
%
% CAUTION: You *must* call InitializePsychSound before first invocation of
% PsychPortAudio(), at least on MS-Windows, but possibly also on OS/X! If
% you omit that call, initialization of the driver may fail with some
% "Invalid MEX file" error from Matlab!
%
%
% PsychPortAudio is built around a modified version of the free, open-source
% PortAudio sound library for portable realtime sound: http://www.portaudio.com

% History:
% 06/07/2007 Written (MK).
% 11/01/2008 Remove warning messages about "early beta release" (MK).

% Some check for not yet supported operating systems:
AssertMex('PsychPortAudio.m');
