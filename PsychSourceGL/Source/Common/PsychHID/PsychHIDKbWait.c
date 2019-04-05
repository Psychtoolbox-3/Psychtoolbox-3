/*
 *   PsychtoolboxGL/Source/Common/PsychHID/PsychHIDWait.c
 *
 *   PROJECTS: None. This function is dead.
 */

#include "PsychHID.h"

PsychError PSYCHHIDKbWait(void) 
{
    PsychErrorExitMsg(PsychError_unimplemented, "This function is no longer supported by PsychHID(). Use KbWait instead.");
    return(PsychError_unimplemented);
}
