/* 
    /osxptb/trunk/PsychSourceGL/Source/OSX/Eyelink/EyelinkReceiveFile.c
 
	PROJECTS: Eyelink 
 
	AUTHORS:
    E.Peters@ai.rug.nl				emp
    f.w.cornelissen@med.rug.nl		fwc
 
	PLATFORMS:	Currently only OS X  
 
	HISTORY:
 
	29-05-2001	emp 		created it
	28/06/06	fwc		adapted from alpha version
 
	TARGET LOCATION:
 
 Eyelink.mexmac resides in:
 PsychHardware/EyelinkToolbox/EyelinkBasic
 
 WARNING, Maybe the implementation of this routine works only on the Macintosh
 
 */
 
#include "PsychEyelink.h"
 
 static char useString[] = "[status =] Eyelink('ReceiveFile',['filename'], ['dest'], ['dest_is_path'])";
 
 static char synopsisString[] = 
 " If <src> is omitted, tracker will send last opened data file.\n"
 " If <dest> is omitted, creates local file with source file name.\n"
 " Else, creates file using <dest> as name.  If <dest_is_path> is supplied and non-zero\n" 
 " uses source file name but adds <dest> as directory path.\n"
 " returns: file size if OK, 0 if file transfer was cancelled, negative =  error code";
 
 static char seeAlsoString[] = "";
 
 /*
  ROUTINE: EyelinkReceiveFile
  PURPOSE:
  uses INT32 receive_data_file(char *src, char *dest, INT16 dest_is_path); 
  If <src> = "", tracker will send last opened data file.
  If <dest> is NULL or "", creates local file with source file name.
  Else, creates file using <dest> as name.  If <dest_is_path> != 0 
  uses source file name but adds <dest> as directory path.
  returns: file size if OK, negative =  error code       */
 
 PsychError EyelinkReceiveFile(void)
 {
	 int iStatus = -1;
	 char *src = "";
	 char *dest = "";
	 int dest_is_path = 0;
	 
	 //all sub functions should have these two lines
	 PsychPushHelp(useString, synopsisString, seeAlsoString);
	 if(PsychIsGiveHelp()){
		 PsychGiveHelp();
		 return(PsychError_none);
	 };
	 
	 //check to see if the user supplied superfluous arguments
	 PsychErrorExit(PsychCapNumInputArgs(3));
	 PsychErrorExit(PsychRequireNumInputArgs(0));
	 PsychErrorExit(PsychCapNumOutputArgs(1));
	 
	 // Verify eyelink is up and running
	 EyelinkSystemIsConnected();
	 EyelinkSystemIsInitialized();
	 
	 
	 PsychAllocInCharArg(1, FALSE, &src);
	 PsychAllocInCharArg(2, FALSE, &dest);
	 PsychCopyInIntegerArg(3, FALSE, &dest_is_path);
	 
	 iStatus = (int) receive_data_file(src, dest, (INT16) dest_is_path);
	 
	 /* if there is an output variable available, assign result to it.   */			
	 PsychCopyOutDoubleArg(1, FALSE, iStatus);
	 
	 return(PsychError_none);
 }
