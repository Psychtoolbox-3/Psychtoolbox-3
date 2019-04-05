/*
    PsychToolbox2/Source/Common/PsychStructGlue.c

    AUTHORS:

    Allen.Ingling@nyu.edu           awi
    mario.kleiner.de@gmail.com      mk

    PLATFORMS: All

    PROJECTS: All

    HISTORY:
    12/31/02  awi   wrote it.
    03/28/11   mk   Make 64-bit clean.
*/

//begin include once 
#ifndef PSYCH_IS_INCLUDED_StructGlue
#define PSYCH_IS_INCLUDED_StructGlue

#include "Psych.h"

psych_bool PsychAllocOutStructArray(    int position, 
                                        PsychArgRequirementType isRequired, 
                                        int numElements,
                                        int numFields, 
                                        const char **fieldNames,  
                                        PsychGenericScriptType **pStruct);

void PsychSetStructArrayStringElement(  const char *fieldName,
                                        int index,
                                        char *text,
                                        PsychGenericScriptType *pStruct);

void PsychSetStructArrayDoubleElement(  const char *fieldName,
                                        int index,
                                        double value,
                                        PsychGenericScriptType *pStruct);

void PsychSetStructArrayBooleanElement( const char *fieldName,
                                        int index,
                                        psych_bool state,
                                        PsychGenericScriptType *pStruct);

//You can use PsychSetStructArrayNativeElement instead of PsychSetStructArrayStructElement, the latter checks that pStructInner is really a struct.   
void PsychSetStructArrayStructElement(  const char *fieldName,
                                        int index,
                                        PsychGenericScriptType *pStructInner,
                                        PsychGenericScriptType *pStructOuter);

void PsychSetStructArrayNativeElement(  const char *fieldName,
                                        int index,
                                        PsychGenericScriptType *nativeElement,
                                        PsychGenericScriptType *pStructOuter);

psych_bool PsychAssignOutStructArray(   int position, 
                                        PsychArgRequirementType isRequired,
                                        PsychGenericScriptType *pStruct);

//end include once
#endif
