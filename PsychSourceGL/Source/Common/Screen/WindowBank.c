/*
	PsychToolbox2/Source/Common/WindowBank.cpp		

	AUTHORS:
	
		Allen.Ingling@nyu.edu		awi 

	PLATFORMS: 
	
		Only OS X for now.


	HISTORY:
	
		07/18/02  awi		Wrote it. This is a second attempt.  The first
							attempted to maintain state with static
							variables internal to functions. Much too 
							awkward for prototyping purposes.  
                07/22/05  mk            Windowbank array is allocated and resized dynamically now, so no limit to maximum
                                        number of windows anymore.

	DESCRIPTION:

		WindowBank contains functions for storing records of open 
		onscreen and offscreen windows.  

		See WindowHelpers for convenience functions which access the fields 
		of records defined in WindowBank.h

		TO DO:

		Change the window index returned to the scripting environment to be the address of the window record.     

*/

#include "Screen.h"


/* ------------------------------------------------------------------ 
	Declare static persistant variables local to WindwBank.cpp variables. 
*/

// MK: windowRecordArrayWINBANK is allocated in InitWindowBank(), released in CloseWindowBank(),
// and dynamically resized in PsychCreateWindowArrayRecord(), if necessary to accomodate more windows.
// We resize in chunks of PSYCH_ALLOC_WINDOW_RECORDS_INC window slots to avoid/reduce possible memory-
// fragmentation and make allocation efficient. This approach is better than using linked-lists, because
// it is easier to implement and still allows access to each windowRecord in constant time via array-lookup,
// avoiding increasing access-delays that would happen with long linked lists.
static PsychWindowRecordType **windowRecordArrayWINBANK=NULL;
static PsychScreenRecordType *screenRecordArrayWINBANK[PSYCH_ALLOC_SCREEN_RECORDS];
static int numWindowRecordsWINBANK=0;

// MK: See InitWindowBank() for initial setup of variables below...
static int PSYCH_MAX_WINDOWS=0; // Maximum number of onscreen-/offscreen-windows/textures...	
static int PSYCH_LAST_WINDOW=0; // the highest possible window index			
static int PSYCH_ALLOC_WINDOW_RECORDS=0; //current length of dynamic array allocated to hold window pointers.	

#define PSYCH_ALLOC_WINDOW_RECORDS_INC 4096  // Increment when extending the window bank...

//Local function prototypes
PsychWindowIndexType FindEmptyWindowIndex(void);


//	Window accessor functions for the outside world. 


/*
	PsychFindScreenWindowFromScreenNumber()
	
	Walk down the list of onscreen windows looking for the window open on the specified screen and and set *winRec to point to its window recored if we find one.
	Otherwise set *winRec to NULL. If screenNumber is kPsychUnaffiliatedWindow, we return the first onscreen window we find, regardless of screen.
*/

void PsychFindScreenWindowFromScreenNumber(int screenNumber, PsychWindowRecordType **winRec)
{
	int							i, numWindows;
	PsychWindowRecordType		**windowArray;
	
	*winRec=NULL;
	PsychCreateVolatileWindowRecordPointerList(&numWindows, &windowArray);
	for(i=0;i<numWindows;i++){
		if(PsychIsOnscreenWindow(windowArray[i])){
			if ((windowArray[i]->screenNumber == screenNumber) || (screenNumber == kPsychUnaffiliatedWindow)) {
				*winRec=windowArray[i];
				break;
			}
		}
	}
	PsychDestroyVolatileWindowRecordPointerList(windowArray);	
} 



/* 
	Initialize list of window record pointers.
*/
PsychError InitWindowBank(void)
{
	PsychWindowIndexType i;
        
        PSYCH_ALLOC_WINDOW_RECORDS=PSYCH_ALLOC_WINDOW_RECORDS_INC; // Initial length of array allocated to hold window pointers.	
    
        // MK: Allocate an initial windowbank of default size PSYCH_ALLOC_WINDOW_RECORDS:
        windowRecordArrayWINBANK=malloc(PSYCH_ALLOC_WINDOW_RECORDS * sizeof(PsychWindowRecordType*));
	if (windowRecordArrayWINBANK==NULL) {
            // Out of memory!
            return(PsychError_outofMemory);
        }
        
        PSYCH_MAX_WINDOWS=PSYCH_ALLOC_WINDOW_RECORDS-10;
        PSYCH_LAST_WINDOW=PSYCH_ALLOC_WINDOW_RECORDS-1; //the highest possible window index			
        
        // Initialize with NULL-Ptrs:
        for(i=PSYCH_FIRST_WINDOW;i<=PSYCH_LAST_WINDOW;i++)
		windowRecordArrayWINBANK[i] = NULL;
        
        return(PsychError_none); //no error
}

 
/*
	Free storage for window records.  All windows should be
	close first. 
*/
PsychError CloseWindowBank(void)
{
	PsychWindowIndexType i;
	
	for(i=PSYCH_FIRST_WINDOW;i<=PSYCH_LAST_WINDOW;i++){
            if(windowRecordArrayWINBANK[i] != NULL)
                    free(windowRecordArrayWINBANK[i]);
        }
	
        // MK: Release the array of pointers itself:
        free(windowRecordArrayWINBANK);

        windowRecordArrayWINBANK=NULL;
        PSYCH_MAX_WINDOWS=0;
        PSYCH_LAST_WINDOW=0;			
        PSYCH_ALLOC_WINDOW_RECORDS=0;	
        
        return(PsychError_none);
}

/* 
	Count the number of windows of the specified type.
*/
int PsychCountOpenWindows(PsychWindowType winType)
{
	int i, total=0;
	
	for(i=PSYCH_FIRST_WINDOW;i<=PSYCH_LAST_WINDOW;i++){
            if(windowRecordArrayWINBANK[i]){
                if(winType==kPsychAnyWindow)
                    ++total;
                else if(windowRecordArrayWINBANK[i]->windowType==winType)
                    ++total;
		}
	}
	return(total);
}


/*
        PsychIsLastWindow()
        
        Is this the last open onscreen window for this screen ?
        SCREENClose needs to know this so that it can release the screen 
*/
psych_bool PsychIsLastOnscreenWindow(PsychWindowRecordType *windowRecord)
{
    int i;

    if(!PsychIsOnscreenWindow(windowRecord))
        return(FALSE);
    for(i=PSYCH_FIRST_WINDOW;i<=PSYCH_LAST_WINDOW;i++){
        if(windowRecordArrayWINBANK[i]){
            if(windowRecordArrayWINBANK[i]->screenNumber == windowRecord->screenNumber  &&  
                windowRecordArrayWINBANK[i]->windowIndex != windowRecord->windowIndex	&&
                (windowRecordArrayWINBANK[i]->windowType == kPsychSingleBufferOnscreen || 
                windowRecordArrayWINBANK[i]->windowType == kPsychDoubleBufferOnscreen))
                return(FALSE);
        }
    }
    return(TRUE);
}


/* 
	Allocate and return an empty window record 
*/
void PsychCreateWindowRecord(PsychWindowRecordType **winRec)
{
        PsychWindowRecordType **tmpwindowRecordArrayWINBANK=NULL;
        PsychWindowIndexType i;
    
        //check for space
        if(numWindowRecordsWINBANK==PSYCH_MAX_WINDOWS) {
            // Windowbank - array is full! We reallocate it, extending it
            // by PSYCH_ALLOC_WINDOW_RECORDS_INC additional slots for additional windows.
            tmpwindowRecordArrayWINBANK=realloc(windowRecordArrayWINBANK, (PSYCH_ALLOC_WINDOW_RECORDS + PSYCH_ALLOC_WINDOW_RECORDS_INC) * sizeof(PsychWindowRecordType*));
            if (tmpwindowRecordArrayWINBANK==NULL) {
                // realloc() failed due to out-of-memory!
                PsychErrorExit(PsychError_outofMemory);   //out of memory
            }

            // Success! Update limits and initialize new slots:
            windowRecordArrayWINBANK = tmpwindowRecordArrayWINBANK;
            PSYCH_ALLOC_WINDOW_RECORDS+=PSYCH_ALLOC_WINDOW_RECORDS_INC;
            PSYCH_MAX_WINDOWS+=PSYCH_ALLOC_WINDOW_RECORDS_INC;
            // Initialize new slots with NULL-Ptrs:
            i=PSYCH_LAST_WINDOW + 1;
            PSYCH_LAST_WINDOW+=PSYCH_ALLOC_WINDOW_RECORDS_INC;
            for(;i<=PSYCH_LAST_WINDOW;i++) windowRecordArrayWINBANK[i] = NULL;
            // Ready for addition of new windows.
        }
    	
	// Allocate storage, zero-init it:
	if((*winRec = (PsychWindowRecordType *) calloc(1, sizeof(PsychWindowRecordType))) == NULL)
		PsychErrorExit(PsychError_outofMemory);   //out of memory
	
	//increment counts	
	++numWindowRecordsWINBANK; 
		
	//store the record at a free pointer index and set the records field to the index.  
	(*winRec)->windowIndex = FindEmptyWindowIndex();
	windowRecordArrayWINBANK[(*winRec)->windowIndex] = *winRec;
        
	//set a flag to indicate that the contents of the window record are not completely valid.
	(*winRec)->isValid=FALSE;
        
	// Set window type to kPsychNoWindow to signal that the window doesn't yet exist.
	// Although the windowRecord is partially initialized, there's no Windowsystem onscreen
	// window associated with this and the OpenGL context isn't there either. This is important
	// for error-handling. Windows of type kPsychNoWindow are ignored by the OpenGL and Window system
	// cleanup routine PsychCloseWindow()...
	(*winRec)->windowType = kPsychNoWindow;

	// Assign default number of color channels: 4 is a good number (RGBA), but this
	// gets overwritten in appropriate places...
	(*winRec)->nrchannels=4;
	
	// Initialize all imaging pipeline related fields to safe defaults:
	PsychInitImagingPipelineDefaultsForWindowRecord(*winRec);
	
	//Intialize the text settings field within the window record to default values which should be (but are not yet) specified in Psychtoolbox preferences.
	PsychInitTextRecordSettings(&((*winRec)->textAttributes));
	
	//Initialize the fields about textures.  
	PsychInitWindowRecordTextureFields(*winRec);
	
	//Initialize the fields used to store alpha blending factors as set by glBlendFunc()
	PsychInitWindowRecordAlphaBlendingFactors(*winRec);
        
	//Initialize line stipple values
	(*winRec)->stipplePattern=0xAAAA;		//alternating pixels stipple pattern
	(*winRec)->stippleFactor=1;
	(*winRec)->stippleEnabled=FALSE;
	
	// Set GL color buffer writemask glColorMask to "all enabled":
	(*winRec)->colorMask[0] = GL_TRUE;
	(*winRec)->colorMask[1] = GL_TRUE;
	(*winRec)->colorMask[2] = GL_TRUE;
	(*winRec)->colorMask[3] = GL_TRUE;
	
	// Initialize stereo settings:
	(*winRec)->stereomode=0;
	(*winRec)->stereodrawbuffer=2;                  // No stero drawbuffer selected at window open time.
	(*winRec)->slaveWindow=NULL;					// No slave window attached.
	(*winRec)->parentWindow=NULL;					// No parent window attached.
	(*winRec)->targetFlipFieldType=-1;				// Don't care if flipping should only happen in even or odd video refresh frames.
	(*winRec)->auxbuffer_dirty[0]=FALSE;            // AUX-Buffers clean on startup.
	(*winRec)->auxbuffer_dirty[1]=FALSE;
	
	// Initialize gamma table pointers for deferred loading:
	(*winRec)->inRedTable = NULL;
	(*winRec)->inGreenTable = NULL;
	(*winRec)->inBlueTable = NULL;
	(*winRec)->loadGammaTableOnNextFlip = 0;
	
	// Set cached display list handles for drawing functions to "uninitialized":
	(*winRec)->fillOvalDisplayList = 0;
	(*winRec)->frameOvalDisplayList = 0;

	// No special flags set by default:
	(*winRec)->specialflags = 0;
	// No capabilities setup yet:
	(*winRec)->gfxcaps = 0;
	
	// No aux param vector set by default:
	(*winRec)->auxShaderParams = NULL;
	(*winRec)->auxShaderParamsCount = 0;

	// Reset memory accounting info for this windowRecord:
	(*winRec)->surfaceSizeBytes = 0;
	
	// NULL out flipinfo struct:
	(*winRec)->flipInfo = NULL;
	
	// Init our shader handles to zero -- Off by default:
	(*winRec)->unclampedDrawShader = 0;
	(*winRec)->defaultDrawShader = 0;

	// Set surface addresses to zero:
	(*winRec)->gpu_preflip_Surfaces[0] = 0;
	(*winRec)->gpu_preflip_Surfaces[1] = 0;

	// Bits per color component value: Initialize to 8 bits per channel. This is true for all
	// standard textures, offscreen- and onscreen windows. Only FBO backed drawables can start
	// off with 16 bpc (float/fixed) or 32 bpc float, or get upgraded to such resolutions. Such
	// upgrades always happen in a maketexture function or FBO creation function.
	// This will be later replaced by true bitdepths for onscreen windows, e.g., to account for
	// 10 bit native framebuffers on latest GPUs:
	(*winRec)->bpc = 8;
	
	// Sync of doublebuffer-swaps to vertical retrace signal starts off as disabled:
	(*winRec)->vSynced = FALSE;
	
	// Default to unknown maximum supported texture size:
	(*winRec)->maxTextureSize = 0;

	// No GPU core id assigned yet:
	(*winRec)->gpuCoreId[0] = 0;
	(*winRec)->gpuRenderTimeQuery = 0;
	(*winRec)->gpuRenderTime = 0.0;

	// No swap group or barrier assigned:
	(*winRec)->swapGroup = 0;
	(*winRec)->swapBarrier = 0;

	// Zero-init OpenML timing values:
	(*winRec)->reference_ust = 0;
	(*winRec)->reference_msc = 0;
	(*winRec)->reference_sbc = 0;
	(*winRec)->osbuiltin_swaptime = 0.0;
	(*winRec)->target_sbc = 0;
	(*winRec)->lastSwaptarget_msc = 0;

	// ...and other timing values:
	(*winRec)->postflip_vbltimestamp = -1.0;

    // Init to zero default:
    (*winRec)->VBL_Startline = 0;
    (*winRec)->VBL_Endline = 0;
    
    // One-Time setup of some parameters needed for userspace OpenGL context at first real bind:
    (*winRec)->needsViewportSetup = TRUE;
    
	return;
}

/*
	Free the window record memory and set the corresponding window array entry to null.
	Its up to the calling function to close the window FIRST.
*/
PsychError FreeWindowRecordFromIndex(PsychWindowIndexType windex)
{	
	
	
	if(windex < PSYCH_FIRST_SCREEN)
		return(PsychError_scumberNotWindex); //I was passed a screen number, not a window index
	if(windex <= PSYCH_LAST_SCREEN)
		return(PsychError_scumberNotWindex); //I was passed a screen number, not a window pointer
	if(windowRecordArrayWINBANK[windex] ==NULL)
		return(PsychError_invalidWindex);    //window does not exist
		
	// Release temporary gamma tables, if any:
	free(windowRecordArrayWINBANK[windex]->inRedTable);
	free(windowRecordArrayWINBANK[windex]->inGreenTable);
	free(windowRecordArrayWINBANK[windex]->inBlueTable);

	free(windowRecordArrayWINBANK[windex]);
	windowRecordArrayWINBANK[windex] = NULL;
	--numWindowRecordsWINBANK;
	return(PsychError_none);
}

/*
	Free the window record memory and set the corresponding window array entry to null.
	Its up to the calling function to close the window FIRST.
*/
PsychError FreeWindowRecordFromPntr(PsychWindowRecordType *winRec)
{
	PsychWindowIndexType windowIndex;
	
	if(winRec == NULL)
		return(PsychError_nullWinRecPntr); //passed a null window pointer
	
	windowIndex = winRec->windowIndex;
	return(FreeWindowRecordFromIndex(windowIndex));
}



/*
	Accept a window pointer or a screen number (psychIndex) and return true if its a valid screen
	number and false otherwise
*/
psych_bool IsValidScreenNumber(PsychNumdexType numdex)
{
	return((int)((numdex >=PSYCH_FIRST_SCREEN) && (numdex < PsychGetNumDisplays())));

}


/*
	IsValidScreenNumberOrUnaffiliated()
	
	Accept a window pointer or a screen number (psychIndex) and return true iff its a valid screen
	number or is the token indicating an unaffiliated window. 
*/
psych_bool IsValidScreenNumberOrUnaffiliated(PsychNumdexType numdex)
{
	return(((int)((numdex >=PSYCH_FIRST_SCREEN) && (numdex < PsychGetNumDisplays()))) || numdex == kPsychUnaffiliatedWindow);

}


/*
	IsScreenNumberUnaffiliated()

*/
psych_bool IsScreenNumberUnaffiliated(PsychNumdexType numdex)
{
	return(numdex == kPsychUnaffiliatedWindow);

}



/*
	Accept a window pointer or a screen number and return true if its a valid window pointer
	number and false otherwise
*/
psych_bool IsWindowIndex(PsychNumdexType numdex)
{


	if(numdex<PSYCH_FIRST_WINDOW || numdex>PSYCH_LAST_WINDOW)
		return(FALSE);  //outside of index range
		
	//check to see if their is window record there.
	if(windowRecordArrayWINBANK[numdex] != NULL)
		return(1);
	
	return(0);		
}


/*
	Accept a window index and return a pointer to the window record, which will be NULL if no
	window exists.  We return an error if the index is invalid.
        
        This function also checks the valid flag in the window to make sure that that fields were
        successfully filled when the window was created. 
        
        All retrieval of window records should go through this function because it checks for valid
        windows.
        
*/
PsychError FindWindowRecord(PsychWindowIndexType windowIndex, PsychWindowRecordType **windowRecord)
{
	// Check for valid index: Must be within bounds of our array of windowRecord pointers, and the referenced slot must be a non-NULL ptr to a windowRecord:
        // Preinit to safe default of NULL:
        *windowRecord = NULL;
	if(windowIndex<PSYCH_FIRST_WINDOW || windowIndex>PSYCH_LAST_WINDOW || ((*windowRecord = windowRecordArrayWINBANK[windowIndex]) == NULL)) return(PsychError_invalidWindex); // Invalid index!

	// It is a windowRecord: Check if it is valid, ie., has been properly initialized by PTB:
	PsychCheckIfWindowRecordIsValid(*windowRecord);  // This would, e.g., fail because of an early exit when the window was created. It is considered an internal error if this triggers.

	// Success.
	return(PsychError_none);
}


/*
	Accept a screen number and return the corresponding screen record
*/
PsychError FindScreenRecord(int screenNumber, PsychScreenRecordType **screenRecord)
{
    if(IsValidScreenNumber(screenNumber)){
            *screenRecord = screenRecordArrayWINBANK[screenNumber];
            return(PsychError_none);
    }else
            return(PsychError_invalidScumber); 	
    
}

/*
    PsychCreateVolatileWindowRecordPointerList()
    
    Allocates memory for and returns an array holding indices of all open windows.  
    
    This is kind of lame.  It iterates over the window list twice.  The window bank stuff needs to be redone anyway.
	
	We don't really have to worry about deallocating this memory because MATLAB will garbage collect it  when 
	the Psychtoolbox call returns.  
    
*/
void PsychCreateVolatileWindowRecordPointerList(int *numWindows, PsychWindowRecordType ***pointerList)
{
    int 			i,j=0;
    PsychWindowRecordType	**tempList; 
    
    *numWindows=0;
    for(i=PSYCH_FIRST_WINDOW;i<=PSYCH_LAST_WINDOW;i++){
        if(windowRecordArrayWINBANK[i])
            ++(*numWindows);
    }
    tempList=(PsychWindowRecordType **)mxMalloc(sizeof(PsychWindowRecordType *) * *numWindows);
    for(i=PSYCH_FIRST_WINDOW;i<=PSYCH_LAST_WINDOW;i++){
        if(windowRecordArrayWINBANK[i])
            tempList[j++]=windowRecordArrayWINBANK[i];
    }
    *pointerList=tempList;     
}


/*
    PsychDestroyVolatileWindowIndexList()
*/
void PsychDestroyVolatileWindowRecordPointerList(PsychWindowRecordType **pointerList)
{
    mxFree((void *)pointerList);
}

/* PsychAssignParentWindow()
 * Assign a specific windowRecord as the parent windowRecord of a windowRecord.
 * This sets up the "parent" pointer and copies a variety of state and settings from the
 * parent to the child, e.g., OpenGL context, device context, color clamp settings etc...,
 * i.e., it performs the one time setup of the childs properties with the parents properties.
 *
 * Usually called by all functions that create new offscreen windows or textures to assign
 * the creating onscreen window to a child offscreen window or texture or proxy.
 *
 */
void  PsychAssignParentWindow(PsychWindowRecordType *childWin, PsychWindowRecordType *parentWin)
{
	// Assign parent:
	childWin->parentWindow = parentWin;
	
	// Copy some state and settings from parent to child:
	memcpy(&childWin->targetSpecific, &parentWin->targetSpecific, sizeof(parentWin->targetSpecific));

	// Copy default drawing shaders from parent:
	childWin->defaultDrawShader   = parentWin->defaultDrawShader;
	childWin->unclampedDrawShader = parentWin->unclampedDrawShader;
	
	// Copy color range and mode from parent window:
	childWin->colorRange = parentWin->colorRange;
	childWin->applyColorRangeToDoubleInputMakeTexture = parentWin->applyColorRangeToDoubleInputMakeTexture;

	// Copy imaging mode flags from parent:
	childWin->imagingMode = parentWin->imagingMode;

	// Inherit capability bits of parents context:
	childWin->gfxcaps = parentWin->gfxcaps;
	childWin->maxTextureSize = parentWin->maxTextureSize;

    // Type of display backend and graphics api:
    childWin->glApiType = parentWin->glApiType;
    childWin->winsysType = parentWin->winsysType;

	return;
}

/* PsychGetParentWindow()
 * Return the windowRecord* of the parent window of a given window, or the windowRecord itself
 * if the window doesn't have a parent because its a top-level onscreen window, ie., it is its
 * own parent.
 */
PsychWindowRecordType* PsychGetParentWindow(PsychWindowRecordType *windowRecord)
{
	while (windowRecord->parentWindow) windowRecord = windowRecord->parentWindow;
	return(windowRecord);
}

//  ------------------------------------------------------------------
//	Accessor functions for stuff internal to WindowBank.cpp.   
//

PsychWindowIndexType FindEmptyWindowIndex(void)
{
	int i;
	for(i=PSYCH_FIRST_WINDOW;i<=PSYCH_LAST_WINDOW;i++){
		if(windowRecordArrayWINBANK[i]==NULL)
			return(i);
	}
	PsychErrorExitMsg(PsychError_toomanyWin,NULL);
	return(PSYCH_INVALID_WINDEX);
}


