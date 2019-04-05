/*
    PsychToolbox3/Source/Common/Screen/RegisterProject.c

    PROJECTS:

        Screen only.

    AUTHORS:

        Allen.Ingling@nyu.edu           awi
        mario.kleiner.de@gmail.com      mk

    PLATFORMS:

        All.

    HISTORY:

        8/20/02  awi    Created.

    DESCRIPTION:

    The psychtoolbox is seperated into two layers:

    -A lower layer which includes C files and all headers included
    by Psych.h.  All Psychtoolbox projects should include Psych.h

    -An upper layer which includes all C files and headers specific
    to the project, for example only the Screen project includes
    Screen.h, which  includes among other things declarations for
    Screen.h functions.

    That the entry point (e.g. MexFunction) for the module is in the
    lower Psych.h layer breaks the typical inclusion heirarchy where the
    library headers are included in the project headers but not vice
    versa. So with Psychtoolbox project the library must include
    the project headers and to do this cleanly and portably we allow
    the lower layer to contain only one header file from the
    upper layer which is alwyas named RegisterProject.h, which in turn
    includes the the header for the project (e.g. Screen.h).  This
    avoids the messiness of having to name project-specific headers
    for every project within a portable library; The library would
    have to be modified to be used in a project.

    Every psychtoolbox project must contain the RegisterProject.h
    header file, and that file should #include the project header file,
    (e.g. Screen.h).

    Every psychtoolbox project should include the function
    PsychModuleInit() which is called when the compiled module is loaded:
    PsychModuleInit should do whatever initilization is the project requires
    and register subfunctions.  The Psychtoolbox provides two functions
    for registering projects:

        PsychError PsychRegisterProjectExit(PsychExitFunctionType *exitFunc)

        -Registers a cleanup function which is invoked before the module
        is purged.

        PsychError PsychRegisterProjectFunction(char *name,  func)

        -When called with the name argument and NULL for the function
        argument it registers name as the project name.

        -When called with name and a non-null func value registers
        func as to be invoked if the module is called with the name as
        the first or second argument.

        -when called with NULL for name and a non-null for func
        registers a module function to be called when the module is
        invoked with no subcommand name.  The resulting behavior varies
        depeinding on whether the project uses named subfunctions.
        The psychoolbox operates in either of two ways depending on
        whether named subfunctions are used.

        -If ONLY the null name function is regsitered then the library
        operates in "direct mode" where control passes to the function
        passed in the func argument whenever the module is invoked.

        -If any named subfunction is registered then the function
        registered with the NULL name is invoked whenever the module
        is invoked with no arguments.

*/

// begin include once

#include "RegisterProject.h"

// This special sauce only on the Unices:
#if PSYCH_SYSTEM != PSYCH_WINDOWS
    // Note: On Linux, _GNU_SOURCE must be defined at the compiler command line, so we
    // can use GNU specific extensions to dlopen(), as well as the dladdr()
    // function below. On OSX, that #define is not needed, as it is part of standard BSD.
    #include <dlfcn.h>
#endif

// PsychModuleInit is in Screen.cpp
// PsychProjectExit is in Screen.cpp
// PsychProjectSelectFunction is in Screen.cpp

PsychError PsychModuleInit(void)
{
    // Initialize all Screen('Preference', ...); settings to startup defaults:
    // This will also set an override level of verbosity if the relevant environment
    // variable has been set for low-level debugging purposes.
    PrepareScreenPreferences();

    // This little Screen reloading gem disabled, as it is bad, bad, bad!
    // Wasn't useful for exporting ConsoleInputHelper(), as that approach proved
    // way too fragile. Additionally it now causes trouble with some workaround for
    // a Mesa bug: Exporting Screen()'s symbols causes moglcore to link against them
    // for GLEW stuff, and thereby it prevents 'clear Screen' if moglcore itself is
    // not cleared itself, as moglcore holds references to symbols in Screen.
    // Disable this code, so symbols of Screen() stay private to Screen.mex itself
    // on Linux and this specific problem is hopefully avoided.
    #if PSYCH_SYSTEM != PSYCH_WINDOWS
    #if 0
        // Try to reopen ourselves, ie. the currently used Screen() MEX file, in "noload" mode,
        // aka RTLD_NOLOAD.
        //
        // We are not really reloaded, as we are obviously already present, but our symbols are
        // made available globally in the process due to reopening with RTLD_GLOBAL, therefore
        // other PTB mex file modules have some access to our functions, allowing cross-module
        // C-Function calls. PsychHID needs this to call our ConsoleInputHelper() function inside
        // SCREENGetMouseHelper.c to reenable keystroke dispatch to octave or matlab -nojvm from
        // their controlling tty's -- To undo a ListenChar(2) issued in such a setup when CTRL+C
        // is detected. There may be future other uses, e.g., with IOPort et al. to trigger or wait
        // completion of actions in Screen() in response to I/O trigger reception, or before emission
        // of I/O trigger/control events.
        //
        // This functionality is Unix specifc, as it is only supported by GLibc with the
        // _GNU_SOURCE flag defined, or on BSD systems like OSX:
        Dl_info screen_info;
        void* screenHandle = NULL;

        // Retrieve the full path and filename of the shared library which contains
        // the C-Function "ScreenExitFunction". As Screen() is the only mex file which
        // defines this function, ie., ourselves, this will hopefully retrieve info
        // about our running loaded instance of Screen():
        if (dladdr((void*) ScreenExitFunction, &screen_info)) {
            // Yes, got the info.
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Screen() plugin path as detected by dladdr(): %s\n", screen_info.dli_fname);

            // Try to re-open us and get a handle to ourselves. This will not really
            // reload us due to RTLD_NOLOAD, just return the handle, increment our
            // refcount and most importantly, make all our public symbols available
            // globally to the whole address space of the runtime process, so other
            // modules can access our non-private symbols and call our non-private
            // functions (aka RTLD_GLOBAL):
            dlerror();
            screenHandle = dlopen(screen_info.dli_fname, RTLD_NOW | RTLD_NOLOAD | RTLD_GLOBAL);
            if (screenHandle) {
                // Worked! Close us again, so the refcount drops back to the old value
                // and the runtime isn't prevented from clear'ing us if needed:
                dlclose(screenHandle);

                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Reopened Screen() mex plugin with global symbol scope!\n");
            }
            else {
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-WARNING: Could not re-dlopen() Screen() from file %s! [%s]\n", screen_info.dli_fname, dlerror());
            }
        } else {
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-WARNING: Could not get module info for Screen() from dladd()!\n");
        }

        if (!screenHandle && (PsychPrefStateGet_Verbosity() > 4)) printf("PTB-WARNING: Other PTB modules can't cross-call Screen's internal functions.\n");
    #endif
    #endif

    // Register the project exit function
    PsychErrorExitMsg(PsychRegisterExit(&ScreenExitFunction), "Failed to register the Screen exit function.");

    // Register the project function which is called when the module
    // is invoked with no arguments:
    PsychErrorExitMsg(PsychRegister(NULL,  &PsychDisplayScreenSynopsis), "Failed to register the Screen synopsis function.");

    // Register the module name
    PsychErrorExitMsg(PsychRegister("Screen", NULL), "Failed to register the Screen module name.");

    // Register named subfunctions

    // Enable for debugging purposes
    PsychErrorExit(PsychRegister("Null",  &SCREENNull));

    // Public api functions:
    PsychErrorExit(PsychRegister("Version",  &MODULEVersion));
    PsychErrorExit(PsychRegister("Computer",  &SCREENComputer));
    PsychErrorExit(PsychRegister("Screens", &SCREENScreens));
    PsychErrorExit(PsychRegister("PixelSize",&SCREENPixelSize));
    PsychErrorExit(PsychRegister("PixelSizes",&SCREENPixelSizes));
    PsychErrorExit(PsychRegister("OpenWindow",  &SCREENOpenWindow));
    PsychErrorExit(PsychRegister("OpenOffscreenWindow",  &SCREENOpenOffscreenWindow));
    PsychErrorExit(PsychRegister("Close",  &SCREENClose));
    PsychErrorExit(PsychRegister("CloseAll",  &SCREENCloseAll));
    PsychErrorExit(PsychRegister("Flip", &SCREENFlip));
    PsychErrorExit(PsychRegister("AsyncFlipBegin", &SCREENFlip));
    PsychErrorExit(PsychRegister("AsyncFlipEnd", &SCREENFlip));
    PsychErrorExit(PsychRegister("AsyncFlipCheckEnd", &SCREENFlip));
    PsychErrorExit(PsychRegister("WaitUntilAsyncFlipCertain" , &SCREENWaitUntilAsyncFlipCertain));
    PsychErrorExit(PsychRegister("FillRect", &SCREENFillRect));
    PsychErrorExit(PsychRegister("GetImage", &SCREENGetImage));
    PsychErrorExit(PsychRegister("PutImage", &SCREENPutImage));
    PsychErrorExit(PsychRegister("HideCursorHelper", &SCREENHideCursorHelper));
    PsychErrorExit(PsychRegister("ShowCursorHelper", &SCREENShowCursorHelper));
    PsychErrorExit(PsychRegister("SetMouseHelper", &SCREENSetMouseHelper));
    PsychErrorExit(PsychRegister("Rect", &SCREENRect));
    PsychErrorExit(PsychRegister("WindowScreenNumber", &SCREENWindowScreenNumber));
    PsychErrorExit(PsychRegister("Windows", &SCREENWindows));
    PsychErrorExit(PsychRegister("WindowKind", &SCREENWindowKind));
    PsychErrorExit(PsychRegister("IsOffscreen", &SCREENIsOffscreen));
    PsychErrorExit(PsychRegister("ReadNormalizedGammaTable", &SCREENReadNormalizedGammaTable));
    PsychErrorExit(PsychRegister("LoadNormalizedGammaTable", &SCREENLoadNormalizedGammaTable));
    PsychErrorExit(PsychRegister("FrameRate", &SCREENNominalFramerate));
    PsychErrorExit(PsychRegister("NominalFrameRate", &SCREENNominalFramerate));
    PsychErrorExit(PsychRegister("glPoint", &SCREENglPoint));
    PsychErrorExit(PsychRegister("gluDisk", &SCREENgluDisk));
    PsychErrorExit(PsychRegister("FillOval", &SCREENFillOval));
    PsychErrorExit(PsychRegister("FrameOval", &SCREENFrameOval));
    PsychErrorExit(PsychRegister("TextModes", &SCREENTextModes));
    PsychErrorExit(PsychRegister("TextMode", &SCREENTextMode));
    PsychErrorExit(PsychRegister("TextSize", &SCREENTextSize));
    PsychErrorExit(PsychRegister("TextStyle", &SCREENTextStyle));
    PsychErrorExit(PsychRegister("TextFont", &SCREENTextFont));
    PsychErrorExit(PsychRegister("TextBounds", &SCREENTextBounds));
    PsychErrorExit(PsychRegister("DrawText", &SCREENDrawText));
    PsychErrorExit(PsychRegister("TextColor", &SCREENTextColor));
    PsychErrorExit(PsychRegister("Preference", &SCREENPreference));
    PsychErrorExit(PsychRegister("MakeTexture", &SCREENMakeTexture));
    PsychErrorExit(PsychRegister("DrawTexture", &SCREENDrawTexture));
    PsychErrorExit(PsychRegister("FrameRect", &SCREENFrameRect));
    PsychErrorExit(PsychRegister("DrawLine", &SCREENDrawLine));
    PsychErrorExit(PsychRegister("FillPoly", &SCREENFillPoly));
    PsychErrorExit(PsychRegister("FramePoly", &SCREENFramePoly));
    PsychErrorExit(PsychRegister("GlobalRect", &SCREENGlobalRect));
    PsychErrorExit(PsychRegister("DrawDots", &SCREENDrawDots));
    PsychErrorExit(PsychRegister("GetTimeList", &SCREENGetTimeList));
    PsychErrorExit(PsychRegister("ClearTimeList", &SCREENClearTimeList));
    PsychErrorExit(PsychRegister("BlendFunction", &SCREENBlendFunction));
    PsychErrorExit(PsychRegister("WindowSize", &SCREENWindowSize));
    PsychErrorExit(PsychRegister("GetMouseHelper", &SCREENGetMouseHelper));
    PsychErrorExit(PsychRegister("TextBackgroundColor", &SCREENTextBackgroundColor));
    PsychErrorExit(PsychRegister("LineStipple", &SCREENLineStipple));
    PsychErrorExit(PsychRegister("SelectStereoDrawBuffer", &SCREENSelectStereoDrawBuffer));
    PsychErrorExit(PsychRegister("DrawingFinished", &SCREENDrawingFinished));
    PsychErrorExit(PsychRegister("DrawLines", &SCREENDrawLines));
    PsychErrorExit(PsychRegister("GetFlipInterval", &SCREENGetFlipInterval));
    PsychErrorExit(PsychRegister("CloseMovie", &SCREENCloseMovie));
    PsychErrorExit(PsychRegister("OpenMovie", &SCREENOpenMovie));
    PsychErrorExit(PsychRegister("PlayMovie", &SCREENPlayMovie));
    PsychErrorExit(PsychRegister("SetMovieTimeIndex", &SCREENSetMovieTimeIndex));
    PsychErrorExit(PsychRegister("GetMovieTimeIndex", &SCREENGetMovieTimeIndex));
    PsychErrorExit(PsychRegister("GetMovieImage", &SCREENGetMovieImage));
    PsychErrorExit(PsychRegister("glPushMatrix", &SCREENglPushMatrix));
    PsychErrorExit(PsychRegister("glPopMatrix", &SCREENglPopMatrix));
    PsychErrorExit(PsychRegister("glLoadIdentity", &SCREENglLoadIdentity));
    PsychErrorExit(PsychRegister("glTranslate", &SCREENglTranslate));
    PsychErrorExit(PsychRegister("glScale", &SCREENglScale));
    PsychErrorExit(PsychRegister("glRotate", &SCREENglRotate));
    PsychErrorExit(PsychRegister("PreloadTextures", &SCREENPreloadTextures));
    PsychErrorExit(PsychRegister("FillArc", &SCREENFillArc));
    PsychErrorExit(PsychRegister("DrawArc", &SCREENDrawArc));
    PsychErrorExit(PsychRegister("FrameArc", &SCREENFrameArc));
    PsychErrorExit(PsychRegister("CopyWindow", &SCREENCopyWindow));
    PsychErrorExit(PsychRegister("WaitBlanking", &SCREENWaitBlanking));
    PsychErrorExit(PsychRegister("GetOpenGLTexture", &SCREENGetOpenGLTexture));
    PsychErrorExit(PsychRegister("SetOpenGLTexture", &SCREENSetOpenGLTexture));
    PsychErrorExit(PsychRegister("BeginOpenGL", &SCREENBeginOpenGL));
    PsychErrorExit(PsychRegister("EndOpenGL", &SCREENEndOpenGL));
    PsychErrorExit(PsychRegister("OpenVideoCapture", &SCREENOpenVideoCapture));
    PsychErrorExit(PsychRegister("CloseVideoCapture", &SCREENCloseVideoCapture));
    PsychErrorExit(PsychRegister("StartVideoCapture", &SCREENStartVideoCapture));
    PsychErrorExit(PsychRegister("StopVideoCapture", &SCREENStopVideoCapture));
    PsychErrorExit(PsychRegister("GetCapturedImage", &SCREENGetCapturedImage));
    PsychErrorExit(PsychRegister("SetVideoCaptureParameter", &SCREENSetVideoCaptureParameter));
    PsychErrorExit(PsychRegister("VideoCaptureDevices", &SCREENVideoCaptureDevices));
    PsychErrorExit(PsychRegister("LoadCLUT", &SCREENLoadCLUT));
    PsychErrorExit(PsychRegister("DisplaySize", &SCREENDisplaySize));
    PsychErrorExit(PsychRegister("SetOpenGLTextureFromMemPointer", &SCREENSetOpenGLTextureFromMemPointer));
    PsychErrorExit(PsychRegister("ColorRange", &SCREENColorRange));
    PsychErrorExit(PsychRegister("HookFunction", &SCREENHookFunction));
    PsychErrorExit(PsychRegister("OpenProxy", &SCREENOpenProxy));
    PsychErrorExit(PsychRegister("TransformTexture", &SCREENTransformTexture));
    PsychErrorExit(PsychRegister("DrawTextures", &SCREENDrawTextures));
    PsychErrorExit(PsychRegister("GetWindowInfo", &SCREENGetWindowInfo));
    PsychErrorExit(PsychRegister("GetOpenGLDrawMode", &SCREENGetOpenGLDrawMode));
    PsychErrorExit(PsychRegister("Resolutions", &SCREENResolutions));
    PsychErrorExit(PsychRegister("Resolution", &SCREENResolution));
    PsychErrorExit(PsychRegister("ConfigureDisplay", &SCREENConfigureDisplay));
    PsychErrorExit(PsychRegister("CreateMovie", &SCREENCreateMovie));
    PsychErrorExit(PsychRegister("FinalizeMovie", &SCREENFinalizeMovie));
    PsychErrorExit(PsychRegister("AddFrameToMovie", &SCREENGetImage));
    PsychErrorExit(PsychRegister("AddAudioBufferToMovie", &SCREENAddAudioBufferToMovie));
    PsychErrorExit(PsychRegister("GetFlipInfo", &SCREENGetFlipInfo));
    PsychErrorExit(PsychRegister("PanelFitter", &SCREENPanelFitter));
    PsychErrorExit(PsychRegister("TextTransform", &SCREENTextTransform));
    PsychErrorExit(PsychRegister("ConstrainCursor", &SCREENConstrainCursor));

    PsychSetModuleAuthorByInitials("awi");
    PsychSetModuleAuthorByInitials("dhb");
    PsychSetModuleAuthorByInitials("dgp");
    PsychSetModuleAuthorByInitials("kas");
    PsychSetModuleAuthorByInitials("fjc");
    PsychSetModuleAuthorByInitials("mk");
    PsychSetModuleAuthorByInitials("cb");

    InitializeSynopsis();
    if (PsychError_none != InitWindowBank())
        printf("PTB-CRITICAL: ERROR: Out of memory when trying to allocate windowBank! Expect complete failure soon!!\n");

    PsychMovieInit();
    PsychVideoCaptureInit();
    PsychMovieWritingInit();

    // Call wait-routine for 0.1 secs: This to initialize the time glue on MS-Windows,
    // so the first call to a timing function won't delay:
    PsychWaitIntervalSeconds(0.1);

    // Reset the "userspaceGL" flag which tells PTB that userspace GL rendering was active
    // due to Screen('BeginOpenGL') command.
    PsychSetUserspaceGLFlag(FALSE);

    // Call display glue init last, as its error handling could go crazy if triggered
    // before the init routines above got executed:
    InitializePsychDisplayGlue();

    return(PsychError_none);
}

#if PSYCH_LANGUAGE == PSYCH_PYTHON

PyMODINIT_FUNC
initScreen(void)
{
    PsychPythonInit("Screen");
}

#endif
