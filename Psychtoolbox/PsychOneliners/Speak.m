function [ ret ] = Speak(saytext, voice, rate, volume, pitch, language)
% Use speech output to speak a given text.
%
% Usage:
%
% [ ret ] = Speak(text [, voice][, rate][, volume][, pitch][, language]);
%
% The function returns an optional 'ret'urn code 0 on success, non-zero
% on failure to speak the requested text.
%
% 'text' must be a text to speak, either a text string or a cell array
% of text strings to speak separately, cell by cell.
%
% The optional 'voice' parameter allows to select among different system
% voices. It is supported on Linux and Mac OS/X.
%
% The names of the available voices differ across operating systems.
%
% Linux supports, e.g., male1,  male2,  male3,  female1,  female2,
% female3, child_male, child_female.
%
% OS/X: Type "!say -v ?" in Matlab to get a list of supported voices.
%
% The optional 'rate' parameter controls speed of speaking on OS/X and
% Linux. On OS/X it defines the number of words per minute, on Linux a
% value between -100 and +100 defines slower or faster speed.
%
% The optional 'volume' parameter allows control of loudness on Linux:
% Value range is -100 to + 100.
%
% The optional 'pitch' parameter allows control of pitch on Linux:
% Value range is -100 to + 100.
%
% The optional 'language' parameter allows control of the output language
% on Linux. E.g., 'de' would output in german language, 'en' english
% language. The text string must be a valid ISO language code string.
%
% Note: Speak on MS-Windows requires the .NET framework to be installed.
% Note: Speak on Linux requires the spd-say command to be installed. This
% is the case by default, e.g., at least on Ubuntu Linux 12.04 and later.
%
% Examples:
% Say "Hello darling" with standard system voice:
% Speak('Hello darling');
%
% Say same text with voice named "Albert":
% Speak('Hello darling', 'Albert');
%

% History:
% 24.07.09 mk           Written for OS/X.
% 03.10.12 Vishal Shah  Added basic support for MS-Windows.
% 06.10.12 mk           Add extended support for OS/X and Linux.

if nargin < 1
    error('You must provide the text string to speak!');
end

% Make saytext cell array of characters:
if ~isa(saytext,'cell')
    saytext = {saytext};
end

if IsOSX
    cmd = 'say ';

    if nargin >= 2 && ~isempty(voice)
        cmd = [cmd sprintf('-v ''%s'' ', voice)];
    end

    if nargin >= 3 && ~isempty(rate)
        cmd = [cmd sprintf('-r %i ', rate)];
    end

    for k=1:length(saytext)
        % Build command string for speech output and do a system() call:
        ret = system(sprintf('%s ''%s''', cmd, saytext{k}));
    end
end

if IsLinux
    cmd = 'spd-say --wait ';

    if nargin >= 2 && ~isempty(voice)
        cmd = [cmd sprintf('--voice-type ''%s'' ', voice)];
    end

    if nargin >= 3 && ~isempty(rate)
        cmd = [cmd sprintf('--rate %i ', rate)];
    end

    if nargin >= 4 && ~isempty(volume)
        cmd = [cmd sprintf('--volume %i ', volume)];
    end

    if nargin >= 5 && ~isempty(pitch)
        cmd = [cmd sprintf('--pitch %i ', pitch)];
    end

    if nargin >= 6 && ~isempty(language)
        cmd = [cmd sprintf('--language ''%s'' ', language)];
    end

    ret = 0;
    for k=1:length(saytext)
        % Build command string for speech output and do a system() call:
        ret = system(sprintf('%s ''%s''', cmd, saytext{k}));
        if ret
            break;
        end
    end

    if ret
        warning('Speak: You need to install the spd-say function (speech-dispatcher) to use this function on Linux. Skipped.'); %#ok<WNTAG>
    end
end

if IsWin
    try
        % Using
        % Microsoft's TTS Namespace
        % http://msdn.microsoft.com/en-us/library/system.speech.synthesis.ttsengine(v=vs.85).aspx
        % Microsoft's Synthesizer Class
        % http://msdn.microsoft.com/en-us/library/system.speech.synthesis.speechsynthesizer(v=vs.85).aspx

        NET.addAssembly('System.Speech');
        Speaker = System.Speech.Synthesis.SpeechSynthesizer;
        for k=1:length(saytext)
            Speaker.Speak (saytext{k});
        end
        ret=0;
    catch
        warning('Speak: You need to install the .Net framework to use this function on Windows. Skipped.'); %#ok<WNTAG>
        ret=1;
    end
end

return;
