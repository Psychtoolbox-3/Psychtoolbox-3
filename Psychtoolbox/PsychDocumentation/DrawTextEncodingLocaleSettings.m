% DrawTextEncodingLocaleSettings
%
% If you are just concerned with drawing "standard" ASCII text with the 
% Screen('DrawText', window, text); command, i.e. text that only contains
% "standard" latin characters, e.g., abcdef....ABCDEF.... and numbers
% 01234567890 and a few punctuation marks and symbols like .,!?#@$%& etc.
% then there's nothing you'd need to do. Screen('DrawText') and
% DrawFormattedText() will process these just fine.
%
% If you want to draw "native" unicode text strings conforming to
% UTF-16 encoding with characters or code points in the BMP (Basic
% multilingual plane), you can simply convert your string of unicode code
% points into a vector of doubles, ie., pass it as a vector of numbers of
% class double(), e.g., text = [1424, 1425, 1426] for the first three
% characters of the hebrew alphabet.
%
% Direct drawing of unicode text characters outside the BMP is possible in the
% same way on GNU/Linux with its default textrenderer, and on OS/X with
% Screen('Preference','TextRenderer', 2); set, if you've installed the
% neccessary software as described in "help DrawTextPlugin". In that case,
% the unicode code points may have values greater than 65535, ie., they are
% processed internally as UTF-32 / UCS-4 characters - Each single number
% represents exactly one character and the numbers can be codepoints from
% any Unicode plane.
%
% The standard OS/X text renderer and MS-Windows text renderer only support
% UTF-16 / UCS-2 encoding. Therefore, characters outside of the BMP must be
% represented by pairs of codepoints, each with a value between 0 and 65535
% - Two consecutive values in your double vector will represent one
% character from outside the BMP.
%
% You can also draw non-ASCII character strings encoded into some single-byte
% encoding scheme, multibyte encoding scheme or unicode characters that are
% not encoded as doubles but as a sequence of UTF-8 encoded bytes. To do
% so, first you have to tell Screen() the encoding scheme of the character
% string by a call to Screen('Preference','TextEncodingLoale', newloc);
% newloc is the name of the text encoding scheme to use. Example names for
% valid locales could be 'en_US.ISO8859-1' for single-byte ISO8859-1
% "Latin-1" encoding, or 'UTF-8' for UTF-8 encoded unicode text, or 'C' for
% default C language character encoding. Passing an empty string '' will
% reset the setting to the system default of your operating system. The
% function returns the previous setting. On OS/X and Linux you can type
% "locale -a" in a terminal window to list all supported settings. For
% MS-Windows you can find supported codepages under this URL:
% http://msdn.microsoft.com/en-us/library/dd317756%28VS.85%29.aspx
% e.g., setting a locale of '.1252' would select MS-Windows codepage 1252,
% aka 'windows-1252'.
% A setting of 'UTF-8' for UTF-8 encoded unicode should always work. The
% set of other supported encodings is highly operating system dependent.
%
% After you've told Screen the encoding to use, you can pass strings to
% Screen('DrawText'). If you want to make sure that Matlab or Octave don't
% tamper in unexpected or uncontrollable ways with your string, you can
% pass the string as a uint8 vector. Passing the string as uint8 vector is
% recommended in most cases. Different versions of Matlab and Octave store
% and process regular char() character strings in different ways on
% different operating systems which can often lead to unexpected results.
% On some setups you will be able to directly pass regular text strings,
% but this is not guaranteed to be portable across different computer
% systems, Matlab or Octave versions.
%
% If you get unexpected results for given text strings you can call
% Screen('Preference','Verbosity', 10); before a drawtext operation. The
% function will then print the passed string as a sequence of numbers after
% conversion into unicode format. This allows to check what kind of unicode
% characters were generated.
%
%
% Assuming Screen('DrawText') managed to receive your final (unicode) text
% string correctly, you'll also have to make sure that a font is selected
% that supports the characters in your string. E.g., not all fonts support
% hebrew characters. Use Screen('TextFont', window, fontName); to select a
% proper font fontName. How to find such a font is operating system
% specific. On Linux, the text renderer allows to select fonts by
% properties instead of names as well, e.g., it allows to select a font
% that supports a requested language:
% Screen('Textfont', window, '-:lang=he');  would select a font that
% supports hebrew. Screen('Textfont', window, '-:lang=ja');  would select a
% font that supports japanese script. See some more examples in the code of
% DrawHighQualityUnicodeTextDemo.m. The same works on OS/X with text
% renderer type 2.
%
% If you type DrawTextEncodingLocaleSettings to run this script, you may
% see some hebrew characters rendered if your systems fonts and text
% renderers do support this. This code also acts as a sample on how to pass
% UTF-8 unicode strings to 'DrawText'.
%

PsychDefaultSetup(1);

try
    screenid = max(Screen('Screens'));
    w = Screen('Openwindow', screenid, 255);
    Screen('TextSize', w, 48);
    if ~IsWin
        Screen('Preference', 'TextRenderer', 2);
        Screen('TextFont', w, '-:lang=he');
    end

    % Some uint8 vector with UTF-8 encoded hebrew characters:
    y = uint8([215  148  215  169  215  160  215  153  215  157   44   32  215  149  215  160 215  162  215  169  215  148]);

    if ~IsLinux
        oldone = Screen('Preference', 'TextEncodingLocale', 'UTF-8');
    else
        oldone = Screen('Preference', 'TextEncodingLocale', 'en_US.UTF-8');
    end

    Screen('DrawText', w, y, 10, 100);
    Screen('Flip', w);

    KbStrokeWait;
    Screen('CloseAll');
    Screen('Preference', 'TextEncodingLocale', oldone);
catch
    sca;

    if exist('oldone','var')
        Screen('Preference', 'TextEncodingLocale', oldone);
    end

    psychrethrow(psychlasterror);
end
