% DrawTextPlugin -- The plugin-based Screen('DrawText') text renderer.
%
% You may have arrived here because Screen() instructed you to go here
% after Screen failed to load the external text renderer plugin.
%
% On both GNU/Linux and Apple MacOS/X, the functions Screen('DrawText') and
% Screen('TextBounds') can use an external text rendering plugin for
% drawing and handling of high quality text. This allows for advanced text
% layout and formatting, high-quality anti-aliased rendering of text at
% arbitrary text sizes, support for modern fonts like TrueType, and support
% for drawing of the full international Unicode character set.
%
% The plugin-based textrenderer can be selected on OS/X and Linux by
% using the command Screen('Preference', 'TextRenderer', 2). On Linux, the
% plugin-based renderer is used by default without need for this command.
% On OS/X, normally Apple's ATSU text renderer would be used, so this
% plugin renderer needs to be manually selected if desired.
%
% The text renderer plugin implements a high-speed renderer based on a
% combination of multiple free software libraries for text rendering and
% text handling:
%
% * OGLFT (http://oglft.sourceforge.net/) the OpenGL-FreeType library.
% * The FreeType-2 (http://freetype.sourceforge.net/) library.
% * The FontConfig (http://www.fontconfig.org) library.
%
% The FontConfig library is used to find the optimal font and font
% settings, given a specific font specification by your user code, a process
% known as "font matching". FontConfig internally uses the FreeType-2
% library to handle the font files on your system and to gather all needed
% information for the matching process.
%
% After a font and settings have been selected, FreeType-2 is used to load
% the font and convert it into high-quality character glyphs, then the
% OGLFT library is used to convert these glyphs into a format optimized for
% fast drawing with OpenGL. OGLFT also performs caching of glyphs, text
% layout, measurement of text dimensions and bounding boxes and the actual
% drawing of the text.
%
% Our actual plugin coordinates all these operations and communicates with
% Screen().
%
% The source code of the plugin can be found in the Psychtoolbox source
% tree under PsychSourceGL/Cohorts/FTGLTextRenderer/
%
% The plugins themselves - one per operating system - are stored in the
% Psychtoolbox/PsychBasic/PsychPlugins folder of your Psychtoolbox
% installation. This is where Screen() expects to find the plugins for
% dynamic loading.
%
% The plugins require a working installation of the FreeType-2 and
% FontConfig libraries on your operating system, ie., somewhere in the
% system library search path.
%
% These libraries are part of the default installation of any decent
% GNU/Linux distribution, so there ain't any need for manual setup work on
% your part.
%
% In the unlikely case they are missing on your Linux setup, something like
% the following commands (on a Debian compatible system like Ubuntu Linux)
% should do the job:
%
% sudo apt-get install freetype
% sudo apt-get install fontconfig
%
% You may need to search around a bit or use your graphical package manager
% if the names don't match perfectly.
%
% On OS/X, these libraries are part of the X-Server, so if you have an
% X-Server installed, you'll have those libraries available.
%
% Reasons for failure to load the plugin:
%
% 1. One of the libraries mentioned above is not installed or not correctly
%    installed on your system. It is also possible that the installed
%    libraries are of an incompatible version.
%
%    cd([PsychtoolboxRoot '/PsychBasic/PsychPlugins']);
%    Run system('ldd libptbdrawtext_ftgl.so.*') on Linux or 
%    Run system('otool -L libptbdrawtext_ftgl.dylib') on OS/X.
%    
%    This will print out a list of all libraries (and their versions)
%    required by the plugin, and some diagnostic to tell you if matching
%    libraries were found or are missing. This should allow you to find out
%    what libraries are missing or incompatible, and to fix the problem.
%
% 2. The plugin can't be found in your Psychtoolbox/PsychBasic/PsychPlugins
%    folder, either because it is missing, or the path to that folder can't
%    be found. In the former case, run UpdatePsychtoolbox again -
%    something's broken with your Psychtoolbox installation. In the latter
%    case you can try to copy the plugin directly into one of the system
%    library folders, e.g., /usr/lib /lib or /usr/local/lib - This way the
%    system may be able to automatically locate the plugin.
%
% 3. Other: Use Google to search the Psychtoolbox forum and the internet
%    for help, check for info on the Psychtoolbox Wiki and - if everything
%    else fails - contact the Psychtoolbox forum for help. In the latter
%    case, describe exactly your system, what troubleshooting steps you
%    performed and the exact results of these steps and their output.
%
