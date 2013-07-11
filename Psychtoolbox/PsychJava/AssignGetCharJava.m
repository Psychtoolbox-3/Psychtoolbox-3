function gcobject = AssignGetCharJava
% AssignGetCharJava - Helper function for Java based GetChar.
%
% This routine is called during initialization from ListenChar,
% CharAvail, GetChar and FlushEvents. It tries as hard as possible
% instantiate a Java object of the GetCharJava class and returns it
% if successfull. Different versions of Matlab use different versions
% of the Java virtual machine. Matlab can only load versions of
% GetCharJava.class that were compiled with the same Java version as
% the one that Matlab uses. Therefore we provide multiple versions of
% the GetCharJava.class, each compiled against a different Java version.
% We try to instantiate the newest version. If that fails, we try the
% 2nd newest version, and so on. If all versions bundled with PTB failed
% to instantiate, the behaviour depends on the availability of an installed
% Java compiler. If a "javac" compiler is available, as on MacOS-X, we invoke
% it on the source code to built a customized version of GetCharJava for that
% machine. Otherwise we give up and fail.

% History:
% 4.10.2006 Written (MK).
% 26.05.2011 Add automatic invocation of PsychJavaTrouble() as the most
%            likely temporary fix for problems (MK).
% 03.07.2013 Drop Java 1.4.2-9 support for pre 10.4 OSX, pre R2007a.
%            Also some code cleanup (MK).

% Try to instantiate different versions of GetCharJava.class, starting with
% the one built against the most recent Java version, descending on failure:
psychlasterror('reset');
try
    % Java version 1.5 -- MacOS-X 10.4.x and M$-Windows Matlab 7.1 and later.
    gcobject = GetCharJava_1_5_0;
    psychlasterror('reset');
    return;
catch %#ok<*CTCH>
    % No op. Just fall through to next case.
end

psychlasterror('reset');
% Using the bundled classes failed. Try to bind a previously built one:
try
    % GetCharJava built from source on this machine via javac compiler.
    gcobject = GetCharJava;
    psychlasterror('reset');
    return;
catch
    % No op. Just fall through to next case.
end

% Check if this is actually a problem of missing PsychJava/ folder on the
% Matlab static classpath. Try to fix it by invocation of
% PsychJavaTrouble():
fprintf('\nWARNING: Failed to load GetCharJava Java class. Most likely your Matlab\n');
fprintf('WARNING: classpath.txt file is not up to date. I will try to auto-fix this for\n');
fprintf('WARNING: the current session by executing the PsychJavaTrouble() command.\n');
fprintf('WARNING: This may have bad side-effects, like clearing all your scripts variables\n');
fprintf('WARNING: and functions! If it works, please read ''help PsychJavaTrouble''\n');
fprintf('WARNING: for a proper and permanent fix to the problem.\n\n');

psychlasterror('reset');
try
    PsychJavaTrouble;

    fprintf('INFO: Problem may be fixed for this session. Let''s see...\n');
    
    gcobject = GetCharJava;
    psychlasterror('reset');
    fprintf('INFO: Yep! Still, please fix the problem properly asap.\n');
    return;
catch
    % No op. Just fall through to next case.
end

% Using the bundled classes failed. Try to invoke the Java compiler on our source
% code to build our own version:
compilecmd = ['javac ' PsychtoolboxRoot 'PsychJava' filesep 'GetCharJava.java'];

fprintf('Could not create a GetCharJava object. This may be due to any of the\n');
fprintf('reasons mentioned in ''help PsychJavaTrouble''. Another reason could be \n');
fprintf('that i do not have a matching GetCharJava implementation for your version\n');
fprintf('of Matlab. Will try now to compile a customized GetCharJava using the javac compiler.\n');
fprintf('This may or may not work...\n');
fprintf('Compile command is: %s\n', compilecmd);

% Invoke compiler, if any:
system(compilecmd);

% Clear out java cache.
psychlasterror('reset');
try
    clear java;
catch
    fprintf('Warning: clear java failed.\n');
end
psychlasterror('reset');

% Retry binding:
try
    % GetCharJava built from source on this machine via javac compiler.
    gcobject = GetCharJava; %#ok<NASGU>
catch
    % Failed again. We are out of luck.
    fprintf('Failed to built customized GetCharJava! This may be due to lack of\n');
    fprintf('an installed Java SDK on your system. Read ''help PsychJavaTrouble'' \n');
    fprintf('for more info.\n\n');
    error('Failed to built and initialize GetCharJava.');
end
psychlasterror('reset');

% Success!
fprintf('Building a GetCharJava file from source on your machine worked!\n');
fprintf('Please exit and restart Matlab to make the new GetCharJava functional.\n');
error('Please exit and restart Matlab.');
