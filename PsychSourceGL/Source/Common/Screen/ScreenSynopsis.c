/*
    PsychToolbox2/Source/Common/Screen/ScreenSynopsis.cpp

    AUTHORS:

        Allen.Ingling@nyu.edu               awi
        mario.kleiner at tuebingen.mpg.de   mk

    PLATFORMS:

        All.

    PROJECTS:

        08/21/02	awi		Screen on MacOS9

    HISTORY:
                08/21/02        awi             Wrote it.
                10/12/04        awi             Added Synopsis strings.
                04/22/05        mk              More/Updated Synopsis strings (OpenWindow, DrawingFinished, SelectStereoDrawbuffer, Flip, DrawLines)
                5/09/05         mk              Added new command GetFlipInterval, changed Flip's synopsis.
                5/13/05         mk              Added synopsis for new rotationAngle argument of "DrawTexture"
                5/30/05         mk              Added synopsis for "SkipSyncTests" - preference setting.
                7/23/05         mk              Added synopsis for new arguments filterMode and globalAlpha of "DrawTexture"
                9/30/05         mk              Added synopsis for "VisualDebugLevel" - preference setting.

        DESCRIPTION:

        Summarize all of screen.  In Matlab this synopsis is given if Screen is
        called with now arguments.


    TO DO:


        ¥The overall Screen Synopsis is  separate from indivisual function help synopses now.  The identical help strings appear within the Screen
        subfunction files and the overal screen synopsis in this file.  It would be better not to duplicate the synopsis strings between files,
        but that gives up the convenience of using the same variable name for the synopsis string in each subfunction file.  There might be a
        clever way to do both using a macro.

        ¥Otherwise, the cool way to do help would be to create an XML DTD for describing
        sub functions which includes their name, short synopsis string,
        long synopsis string, and function arguments and their types.

        From an XML file we then-
        -Validate arguments at runtime.
        -Issue verbose warnings when passed invalid arguments.
        -Parse arguments at runtime.
        -Automatically generate HTML documentation
        -Generate built-in screen help.

        ¥The error issued in InitializeSynopsis is issued after
        the fatal mistake is made not before. Probably not a big deal: because
        Screen always does this init there is no test of Screen which would fail to detect the problem.

*/


#include "Screen.h"

#define MAX_SYNOPSIS_STRINGS 500  

//declare variables local to this file.  
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

void InitializeSynopsis()
{
    int i=0;
    const char **synopsis = synopsisSYNOPSIS;  //abbreviate the long name

    synopsis[i++] = "Usage:";

    // OS-9 PTB emulation:
    synopsis[i++] = "\n% Activate compatibility mode: Try to behave like the old MacOS-9 Psychtoolbox:";
    synopsis[i++] = "oldEnableFlag=Screen('Preference', 'EmulateOldPTB', [enableFlag]);";

    // Open or close a window or texture:
    synopsis[i++] = "\n% Open or close a window or texture:";
    synopsis[i++] = "[windowPtr,rect]=Screen('OpenWindow',windowPtrOrScreenNumber [,color] [,rect] [,pixelSize] [,numberOfBuffers] [,stereomode] [,multisample][,imagingmode][,specialFlags][,clientRect]);";
    synopsis[i++] = "[windowPtr,rect]=Screen('OpenOffscreenWindow',windowPtrOrScreenNumber [,color] [,rect] [,pixelSize] [,specialFlags] [,multiSample]);";
    synopsis[i++] = "textureIndex=Screen('MakeTexture', WindowIndex, imageMatrix [, optimizeForDrawAngle=0] [, specialFlags=0] [, floatprecision=0] [, textureOrientation=0] [, textureShader=0]);";
    synopsis[i++] = "oldParams = Screen('PanelFitter', windowPtr [, newParams]);";
    synopsis[i++] = "Screen('Close', [windowOrTextureIndex or list of textureIndices/offscreenWindowIndices]);";
    synopsis[i++] = "Screen('CloseAll');";

    // Draw lines lines solids like QuickDraw and DirectX (OS 9 and Windows)
    synopsis[i++] = "\n%  Draw lines and solids like QuickDraw and DirectX (OS 9 and Windows):";
    synopsis[i++] = "currentbuffer = Screen('SelectStereoDrawBuffer', windowPtr [, bufferid] [, param1]);";
    synopsis[i++] = "Screen('DrawLine', windowPtr [,color], fromH, fromV, toH, toV [,penWidth]);";
    synopsis[i++] = "Screen('DrawArc',windowPtr,[color],[rect],startAngle,arcAngle)";
    synopsis[i++] = "Screen('FrameArc',windowPtr,[color],[rect],startAngle,arcAngle[,penWidth] [,penHeight] [,penMode])";
    synopsis[i++] = "Screen('FillArc',windowPtr,[color],[rect],startAngle,arcAngle)";
    synopsis[i++] = "Screen('FillRect', windowPtr [,color] [,rect] );";
    synopsis[i++] = "Screen('FrameRect', windowPtr [,color] [,rect] [,penWidth]);";
    synopsis[i++] = "Screen('FillOval', windowPtr [,color] [,rect] [,perfectUpToMaxDiameter]);";
    synopsis[i++] = "Screen('FrameOval', windowPtr [,color] [,rect] [,penWidth] [,penHeight] [,penMode]);";
    synopsis[i++] = "Screen('FramePoly', windowPtr [,color], pointList [,penWidth]);";
    synopsis[i++] = "Screen('FillPoly', windowPtr [,color], pointList [, isConvex]);";

    // New OpenGL-based functions for OS X
    synopsis[i++] = "\n% New OpenGL functions for OS X:";
    synopsis[i++] = "Screen('glPoint', windowPtr, color, x, y [,size]);";
    synopsis[i++] = "Screen('gluDisk', windowPtr, color, x, y [,size]);";
    synopsis[i++] = "[minSmoothPointSize, maxSmoothPointSize, minAliasedPointSize, maxAliasedPointSize] = Screen('DrawDots', windowPtr, xy [,size] [,color] [,center] [,dot_type] [,lenient]);";
    synopsis[i++] = "[minSmoothLineWidth, maxSmoothLineWidth, minAliasedLineWidth, maxAliasedLineWidth] = Screen('DrawLines', windowPtr, xy [,width] [,colors] [,center] [,smooth] [,lenient]);";
    synopsis[i++] = "[sourceFactorOld, destinationFactorOld, colorMaskOld]=Screen('BlendFunction', windowIndex, [sourceFactorNew], [destinationFactorNew], [colorMaskNew]);";

    // Draw Text in windows
    synopsis[i++] = "\n% Draw Text in windows";
    synopsis[i++] = "textModes = Screen('TextModes');";
    synopsis[i++] = "oldCopyMode=Screen('TextMode', windowPtr [,textMode]);";
    synopsis[i++] = "oldTextSize=Screen('TextSize', windowPtr [,textSize]);";
    synopsis[i++] = "oldStyle=Screen('TextStyle', windowPtr [,style]);";
    synopsis[i++] = "[oldFontName,oldFontNumber,oldTextStyle]=Screen('TextFont', windowPtr [,fontNameOrNumber][,textStyle]);";
    synopsis[i++] = "[normBoundsRect, offsetBoundsRect, textHeight, xAdvance] = Screen('TextBounds', windowPtr, text [,x] [,y] [,yPositionIsBaseline] [,swapTextDirection]);";
    synopsis[i++] = "[newX, newY, textHeight]=Screen('DrawText', windowPtr, text [,x] [,y] [,color] [,backgroundColor] [,yPositionIsBaseline] [,swapTextDirection]);";
    synopsis[i++] = "oldTextColor=Screen('TextColor', windowPtr [,colorVector]);";
    synopsis[i++] = "oldTextBackgroundColor=Screen('TextBackgroundColor', windowPtr [,colorVector]);";
    synopsis[i++] = "oldMatrix = Screen('TextTransform', windowPtr [, newMatrix]);";

    // Copy an image, very quickly, between textures and onscreen windows
    synopsis[i++] = "\n% Copy an image, very quickly, between textures, offscreen windows and onscreen windows.";
    synopsis[i++] = "[resident [texidresident]] = Screen('PreloadTextures', windowPtr [, texids]);";
    synopsis[i++] = "Screen('DrawTexture', windowPointer, texturePointer [,sourceRect] [,destinationRect] [,rotationAngle] [, filterMode] [, globalAlpha] [, modulateColor] [, textureShader] [, specialFlags] [, auxParameters]);";
    synopsis[i++] = "Screen('DrawTextures', windowPointer, texturePointer(s) [, sourceRect(s)] [, destinationRect(s)] [, rotationAngle(s)] [, filterMode(s)] [, globalAlpha(s)] [, modulateColor(s)] [, textureShader] [, specialFlags] [, auxParameters]);";
    synopsis[i++] = "Screen('CopyWindow', srcWindowPtr, dstWindowPtr, [srcRect], [dstRect], [copyMode])";

    // Copy an image, slowly, between matrices and windows
    synopsis[i++] = "\n% Copy an image, slowly, between matrices and windows :";
    synopsis[i++] = "imageArray=Screen('GetImage', windowPtr [,rect] [,bufferName] [,floatprecision=0] [,nrchannels=3])";
    synopsis[i++] = "Screen('PutImage', windowPtr, imageArray [,rect]);";

    // Synchronize with the window's screen (on-screen only):
    synopsis[i++] = "\n% Synchronize with the window's screen (on-screen only):";
    synopsis[i++] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', windowPtr [, when] [, dontclear] [, dontsync] [, multiflip]);";
    synopsis[i++] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipBegin', windowPtr [, when] [, dontclear] [, dontsync] [, multiflip]);";
    synopsis[i++] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipEnd', windowPtr);";
    synopsis[i++] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipCheckEnd', windowPtr);";
    synopsis[i++] = "[VBLTimestamp StimulusOnsetTime swapCertainTime] = Screen('WaitUntilAsyncFlipCertain', windowPtr);";
    synopsis[i++] = "[info] = Screen('GetFlipInfo', windowPtr [, infoType=0] [, auxArg1]);";
    synopsis[i++] = "[telapsed] = Screen('DrawingFinished', windowPtr [, dontclear] [, sync]);";
    synopsis[i++] = "framesSinceLastWait = Screen('WaitBlanking', windowPtr [, waitFrames]);";

    // Load color lookup table of the window's screen (on-screen only)
    synopsis[i++] = "\n% Load color lookup table of the window's screen (on-screen only):";
    synopsis[i++] = "[gammatable, dacbits, reallutsize] = Screen('ReadNormalizedGammaTable', windowPtrOrScreenNumber [, physicalDisplay]);";
    synopsis[i++] = "[oldtable, success] = Screen('LoadNormalizedGammaTable', windowPtrOrScreenNumber, table [, loadOnNextFlip][, physicalDisplay][, ignoreErrors]);";
    synopsis[i++] = "oldclut = Screen('LoadCLUT', windowPtrOrScreenNumber [, clut] [, startEntry=0] [, bits=8]);";

    // Get and set information about a window or screen.
    synopsis[i++] = "\n% Get (and set) information about a window or screen:";
    synopsis[i++] = "screenNumbers=Screen('Screens' [, physicalDisplays]);";
    synopsis[i++] = "windowPtrs=Screen('Windows');";
    synopsis[i++] = "kind=Screen(windowPtr, 'WindowKind');";
    synopsis[i++] = "isOffscreen=Screen(windowPtr,'IsOffscreen');";
    synopsis[i++] = "hz=Screen('FrameRate', windowPtrOrScreenNumber [, mode] [, reqFrameRate]);";
    synopsis[i++] = "hz=Screen('NominalFrameRate', windowPtrOrScreenNumber [, mode] [, reqFrameRate]);";
    synopsis[i++] = "[ monitorFlipInterval nrValidSamples stddev ]=Screen('GetFlipInterval', windowPtr [, nrSamples] [, stddev] [, timeout]);";
    synopsis[i++] = "screenNumber=Screen('WindowScreenNumber', windowPtr);";
    synopsis[i++] = "rect=Screen('Rect', windowPtrOrScreenNumber [, realFBSize=0]);";
    synopsis[i++] = "pixelSize=Screen('PixelSize', windowPtrOrScreenNumber);";
    synopsis[i++] = "pixelSizes=Screen('PixelSizes', windowPtrOrScreenNumber);";
    synopsis[i++] = "[width, height]=Screen('WindowSize', windowPointerOrScreenNumber [, realFBSize=0]);";
    synopsis[i++] = "[width, height]=Screen('DisplaySize', ScreenNumber);";
    synopsis[i++] = "[oldmaximumvalue, oldclampcolors, oldapplyToDoubleInputMakeTexture] = Screen('ColorRange', windowPtr [, maximumvalue][, clampcolors][, applyToDoubleInputMakeTexture]);";
    synopsis[i++] = "info = Screen('GetWindowInfo', windowPtr [, infoType=0] [, auxArg1]);";
    synopsis[i++] = "resolutions=Screen('Resolutions', screenNumber);";
    synopsis[i++] = "oldResolution=Screen('Resolution', screenNumber [, newwidth] [, newheight] [, newHz] [, newPixelSize] [, specialMode]);";
    synopsis[i++] = "oldSettings = Screen('ConfigureDisplay', setting, screenNumber, outputId [, newwidth][, newheight][, newHz][, newX][, newY]);";

    // Get and set information about the environment, computer, and video card (i.e. screen):
    synopsis[i++] = "\n% Get/set details of environment, computer, and video card (i.e. screen):";
    synopsis[i++] = "struct=Screen('Version');";
    synopsis[i++] = "comp=Screen('Computer');";
    synopsis[i++] = "oldBool=Screen('Preference', 'IgnoreCase' [,bool]);";
    synopsis[i++] = "tick0Secs=Screen('Preference', 'Tick0Secs', tick0Secs);";
    synopsis[i++] = "psychTableVersion=Screen('Preference', 'PsychTableVersion');";
    synopsis[i++] = "mexFunctionName=Screen('Preference', 'PsychTableCreator');";
    synopsis[i++] = "proc=Screen('Preference', 'Process');";
    synopsis[i++] = "Screen('Preference','SkipSyncTests', skipTest);";
    synopsis[i++] = "Screen('Preference','VisualDebugLevel', level (valid values between 0 and 5));";
    synopsis[i++] = "Screen('Preference', 'ConserveVRAM', mode (valid values between 0 and 3));";
    synopsis[i++] = "Screen('Preference', 'Enable3DGraphics', [enableFlag]);";


    //synopsis[i++] = "\n% Set clipping region (on- or off- screen):";

    // Helper functions.  Don't call these directly, use eponymous wrappers.
    synopsis[i++] = "\n% Helper functions.  Don't call these directly, use eponymous wrappers:";
    synopsis[i++] ="[x, y, buttonVector, hasKbFocus, valuators]= Screen('GetMouseHelper', numButtons [, screenNumber][, mouseIndex]);";
    synopsis[i++] = "Screen('HideCursorHelper', windowPntr [, mouseIndex]);";
    synopsis[i++] = "Screen('ShowCursorHelper', windowPntr [, cursorshapeid][, mouseIndex]);";
    synopsis[i++] = "Screen('SetMouseHelper', windowPntrOrScreenNumber, x, y [, mouseIndex][, detachFromMouse]);";

    // Internal testing of Screen
    synopsis[i++] = "\n% Internal testing of Screen";
    synopsis[i++] =  "timeList= Screen('GetTimelist');";
    synopsis[i++] =  "Screen('ClearTimelist');";
    synopsis[i++] =  "Screen('Preference','DebugMakeTexture', enableDebugging);";

    // Movie and multimedia handling functions:
    synopsis[i++] = "\n% Movie and multimedia playback functions:";
    synopsis[i++] =  "[ moviePtr [duration] [fps] [width] [height] [count] [aspectRatio]]=Screen('OpenMovie', windowPtr, moviefile [, async=0] [, preloadSecs=1] [, specialFlags1=0][, pixelFormat=4][, maxNumberThreads=-1][, movieOptions]);";
    synopsis[i++] =  "Screen('CloseMovie', moviePtr);";
    synopsis[i++] =  "[ texturePtr [timeindex]]=Screen('GetMovieImage', windowPtr, moviePtr, [waitForImage], [fortimeindex], [specialFlags = 0] [, specialFlags2 = 0]);";
    synopsis[i++] =  "[droppedframes] = Screen('PlayMovie', moviePtr, rate, [loop], [soundvolume]);";
    synopsis[i++] =  "timeindex = Screen('GetMovieTimeIndex', moviePtr);";
    synopsis[i++] =  "[oldtimeindex] = Screen('SetMovieTimeIndex', moviePtr, timeindex [, indexIsFrames=0]);";
    synopsis[i++] =  "moviePtr = Screen('CreateMovie', windowPtr, movieFile [, width][, height][, frameRate=30][, movieOptions][, numChannels=4][, bitdepth=8]);";
    synopsis[i++] =  "Screen('FinalizeMovie', moviePtr);";
    synopsis[i++] =  "Screen('AddFrameToMovie', windowPtr [,rect] [,bufferName] [,moviePtr=0] [,frameduration=1]);";
    synopsis[i++] =  "Screen('AddAudioBufferToMovie', moviePtr, audioBuffer);";

    // Video capture support:
    synopsis[i++] = "\n% Video capture functions:";
    synopsis[i++] = "devices = Screen('VideoCaptureDevices' [, engineId]);";
    synopsis[i++] = "videoPtr =Screen('OpenVideoCapture', windowPtr [, deviceIndex][, roirectangle][, pixeldepth][, numbuffers][, allowfallback][, targetmoviename][, recordingflags][, captureEngineType][, bitdepth=8]);";
    synopsis[i++] = "Screen('CloseVideoCapture', capturePtr);";
    synopsis[i++] = "[fps starttime] = Screen('StartVideoCapture', capturePtr [, captureRateFPS] [, dropframes=0] [, startAt]);";
    synopsis[i++] = "droppedframes = Screen('StopVideoCapture', capturePtr [, discardFrames=1]);";
    synopsis[i++] = "[ texturePtr [capturetimestamp] [droppedcount] [average_intensityOrRawImageMatrix]]=Screen('GetCapturedImage', windowPtr, capturePtr [, waitForImage=1] [,oldTexture] [,specialmode] [,targetmemptr]);";
    synopsis[i++] = "oldvalue = Screen('SetVideoCaptureParameter', capturePtr, 'parameterName' [, value]);";

    // Low level OpenGL calls - directly translated to C via very thin wrapper functions:
    synopsis[i++] = "\n% Low level direct access to OpenGL-API functions:";
    synopsis[i++] = "% Online info for each function available by opening a terminal window";
    synopsis[i++] = "% and typing 'man Functionname' + Enter.\n";
    synopsis[i++] = "Screen('glPushMatrix', windowPtr);";
    synopsis[i++] = "Screen('glPopMatrix', windowPtr);";
    synopsis[i++] = "Screen('glLoadIdentity', windowPtr);";
    synopsis[i++] = "Screen('glTranslate', windowPtr, tx, ty [, tz]);";
    synopsis[i++] = "Screen('glScale', windowPtr, sx, sy [, sz]);";
    synopsis[i++] = "Screen('glRotate', windowPtr, angle, [rx = 0], [ry = 0] ,[rz = 1]);";

    // Interfacing with external OpenGL rendering code (MOGL and external OpenGL Mexfiles):
    synopsis[i++] = "\n% Support for 3D graphics rendering and for interfacing with external OpenGL code:";
    synopsis[i++] = "Screen('Preference', 'Enable3DGraphics', [enableFlag]);  % Enable 3D gfx support.";
    synopsis[i++] = "Screen('BeginOpenGL', windowPtr [, sharecontext]);  % Prepare window for external OpenGL drawing.";
    synopsis[i++] = "Screen('EndOpenGL', windowPtr);  % Finish external OpenGL drawing.";
    synopsis[i++] = "[targetwindow, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');";
    synopsis[i++] = "[textureHandle rect] = Screen('SetOpenGLTextureFromMemPointer', windowPtr, textureHandle, imagePtr, width, height, depth [, upsidedown][, target][, glinternalformat][, gltype][, extdataformat][, specialFlags]);";
    synopsis[i++] = "[textureHandle rect] = Screen('SetOpenGLTexture', windowPtr, textureHandle, glTexid, target [, glWidth][, glHeight][, glDepth][, textureShader][, specialFlags]);";
    synopsis[i++] = "[ gltexid gltextarget texcoord_u texcoord_v ] =Screen('GetOpenGLTexture', windowPtr, textureHandle [, x][, y]);";

    synopsis[i++] = "\n% Support for plugins and for builtin high performance image processing pipeline:";
    synopsis[i++] = "[ret1, ret2, ...] = Screen('HookFunction', windowPtr, 'Subcommand', 'HookName', arg1, arg2, ...);";
    synopsis[i++] = "proxyPtr = Screen('OpenProxy', windowPtr [, imagingmode]);";
    synopsis[i++] = "transtexid = Screen('TransformTexture', sourceTexture, transformProxyPtr [, sourceTexture2][, targetTexture][, specialFlags]);";

    synopsis[i++] = NULL;  //this tells PsychDisplayScreenSynopsis where to stop

    if (i > MAX_SYNOPSIS_STRINGS) {
        PrintfExit("%s: increase dimension of synopsis[] from %ld to at least %ld and recompile.",__FILE__,(long)MAX_SYNOPSIS_STRINGS,(long)i);
    }
}

PsychError PsychDisplayScreenSynopsis(void)
{
    int i;

    for (i = 0; synopsisSYNOPSIS[i] != NULL; i++)
        printf("%s\n",synopsisSYNOPSIS[i]);

    return(PsychError_none);
}
