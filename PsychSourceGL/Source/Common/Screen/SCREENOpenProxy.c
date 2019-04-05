/*
    SCREENOpenProxy.c

    AUTHORS:

        mario.kleiner.de@gmail.com      mk

    PLATFORMS:

        All.

    HISTORY:

        01/25/07    mk        Wrote it.

    DESCRIPTION:

        Create a proxy object. A proxy is a windowRecord (like a texture, offscreen window or
        onscreen window), but without associated image content. Its useful as a container for
        all the typical settings associated with a window or texture. Its main purpose is to
        store additional hook function chains for use with the imaging pipeline. The users
        script can create special proxy objects for different image processing subtasks, attach
        proper user defined hook chains to it, and then switch between different hook chains by
        simply passing different proxy objects to the image processing functions. This way we
        are not limited to a fixed number of imaging hooks but can create new ones dynamically.
*/

#include "Screen.h"

// If you change useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[] = "proxyPtr = Screen('OpenProxy', windowPtr [, imagingmode]);";
//                                                        1            2
static char synopsisString[] =
        "Create a proxy object. A proxy is a windowRecord (like a texture, offscreen window or "
        "onscreen window), but without associated image content. Its useful as a container for "
        "all the typical settings associated with a window or texture. Its main purpose is to "
        "store additional hook function chains for use with the imaging pipeline. The users "
        "script can create special proxy objects for different image processing subtasks, attach "
        "proper user defined hook chains to it, and then switch between different hook chains by "
        "simply passing different proxy objects to the image processing functions. This way we "
        "are not limited to a fixed number of imaging hooks but can create new ones dynamically. "
        "Read 'help PsychGLImageProcessing' for more infos about proxies and the pipeline. "
        "\n\n"
        "'windowPtr' is a handle to an associated parent onscreen window. Just pass the handle of "
        "the onscreen window for which image processing is supposed to take place. The routine "
        "returns a handle 'proxyPtr' to the new proxy object. You can delete proxy objects via "
        "the regular Screen('Close', proxyPtr) or Screen('CloseAll'); command as you would do "
        "with textures or windows. You can add or modify hook functions by passing the 'proxyPtr' "
        "to the Screen('Hookfunction', ...) subfunction. You can apply a proxy via the Screen('TransformTexture') "
        "function. 'imagingmode' flags to define the boundary conditions for image processing. ";

static char seeAlsoString[] = "Hookfunction, TransformTexture";

PsychError SCREENOpenProxy(void)
{
    PsychWindowRecordType *windowRecord, *proxyRecord;
    int imagingmode;

    // All subfunctions should have these two lines.
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));
    PsychErrorExit(PsychRequireNumInputArgs(1));
    PsychErrorExit(PsychCapNumOutputArgs(1));

    // Get the window structure for the onscreen window.
    PsychAllocInWindowRecordArg(1, kPsychArgRequired, &windowRecord);

    // Create proxy object:
    PsychCreateWindowRecord(&proxyRecord);  // This also fills the window index field.

    // Set type:
    proxyRecord->windowType=kPsychProxyWindow;

    // Assign parent window and copy its inheritable properties:
    PsychAssignParentWindow(proxyRecord, windowRecord);

    // Init proxy's rects to parent windows. We could set anything that is a non-empty
    // rect, doesn't need to be the windowRecord's rect, but then why not?
    // Without this, functions like PsychSetupView() might fail due to empty
    // rects triggering OpenGL errors:
    PsychCopyRect(proxyRecord->rect, windowRecord->rect);
    PsychCopyRect(proxyRecord->clientrect, windowRecord->clientrect);

    // Get optional imagingmode argument and assign it:
    imagingmode = 0;
    PsychCopyInIntegerArg(2, FALSE, &imagingmode);
    if (imagingmode < 0)
        PsychErrorExitMsg(PsychError_user, "'imagingmode' argument must be a positive number (>=0).");

    proxyRecord->imagingMode = imagingmode;

    // If this will ever get a backing texture, it will be upright/normalized:
    proxyRecord->textureOrientation = 2;

    // Window ready. Mark it valid and return handle to userspace:
    PsychSetWindowRecordValid(proxyRecord);

    //Return the window index and the rect argument.
    PsychCopyOutDoubleArg(1, FALSE, proxyRecord->windowIndex);

    // Done.
    return(PsychError_none);
}
