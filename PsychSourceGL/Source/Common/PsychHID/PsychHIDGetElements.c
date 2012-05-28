/*
  PsychToolbox3/Source/Common/PsychHID/PsychHIDGetElements.c		
  
  PROJECTS: PsychHID
  
  PLATFORMS:  OSX  
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
      
  HISTORY:
  4/29/03  awi		Created.
  
  TO DO:
  

*/

#include "PsychHID.h"

#if PSYCH_SYSTEM == PSYCH_OSX

static char useString[]= "elements=PsychHID('Elements',deviceNumber)";
static char synopsisString[] = 
        "Return a struct array describing each element of a connected USB HID device.";
static char seeAlsoString[] = "";


PsychError PSYCHHIDGetElements(void) 
{
    pRecDevice 			specDevice=NULL;
    UInt32              numDeviceElements;
    
    const char 			*elementFieldNames[]={"typeMaskName", "name", "deviceIndex", "elementIndex", "typeValue", "typeName", "usagePageValue", "usageValue", 							"usageName", "dataSize", "rangeMin", "rangeMax", "scaledRangeMin", "scaledRangeMax", "relative", 
                                                        "wrapping", "nonLinear", "preferredState", "nullState", "calMin", "calMax", "scalingMin", "scalingMax"};
    int 			numElementStructElements, numElementStructFieldNames=23, elementIndex, deviceIndex;
    PsychGenericScriptType	*elementStruct;	
    pRecElement     currentElement, lastElement = NULL;
    char			elementTypeName[PSYCH_HID_MAX_DEVICE_ELEMENT_TYPE_NAME_LENGTH];	
    char			usageName[PSYCH_HID_MAX_DEVICE_ELEMENT_USAGE_NAME_LENGTH];
    char            tmpName[1024];
    char			*typeMaskName;
    HIDElementTypeMask		typeMask;
    	 
    //all subfunctions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
    
    PsychCopyInIntegerArg(1, TRUE, &deviceIndex);
    PsychHIDVerifyInit();
    specDevice= PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
    PsychHIDVerifyOpenDeviceInterfaceFromDeviceIndex(deviceIndex);
    numDeviceElements= HIDCountDeviceElements(specDevice, kHIDElementTypeIO);
    numElementStructElements = (int)numDeviceElements;
    PsychAllocOutStructArray(1, FALSE, numElementStructElements, numElementStructFieldNames, elementFieldNames, &elementStruct);
    elementIndex=0;
    for(currentElement=HIDGetFirstDeviceElement(specDevice,kHIDElementTypeIO); 
        (currentElement != NULL) && (currentElement != lastElement);
        currentElement=HIDGetNextDeviceElement(currentElement, kHIDElementTypeIO)) {
        lastElement = currentElement;

#ifdef __LP64__
        // Skip this entry if it isn't a valid element:
        if (!HIDIsValidElement(currentElement)) continue;

        IOHIDElementType type = IOHIDElementGetType(currentElement);
        typeMask = HIDConvertElementTypeToMask(type);
        PsychHIDGetTypeMaskStringFromTypeMask(typeMask, &typeMaskName);
        PsychSetStructArrayStringElement("typeMaskName",	elementIndex, 	typeMaskName,	 	elementStruct);
        tmpName[0] = 0;
        CFStringRef cfString = IOHIDElementGetName(currentElement);
        if (cfString) {
            CFStringGetCString(cfString, tmpName, sizeof(tmpName), kCFStringEncodingASCII);
            CFRelease(cfString);
        }
        PsychSetStructArrayStringElement("name",            elementIndex, 	tmpName,	 		elementStruct);
        PsychSetStructArrayDoubleElement("deviceIndex",		elementIndex, 	(double)deviceIndex, 			elementStruct);
        PsychSetStructArrayDoubleElement("elementIndex",	elementIndex, 	(double)elementIndex+1, 		elementStruct);
        
        PsychSetStructArrayDoubleElement("typeValue",		elementIndex, 	(double) type, 		elementStruct);
        HIDGetTypeName(type, elementTypeName);
        PsychSetStructArrayStringElement("typeName",		elementIndex, 	elementTypeName,	elementStruct);
        PsychSetStructArrayDoubleElement("usagePageValue",	elementIndex, 	(double)IOHIDElementGetUsagePage(currentElement), 	elementStruct);
        PsychSetStructArrayDoubleElement("usageValue",		elementIndex, 	(double)IOHIDElementGetUsage(currentElement), 		elementStruct);
        HIDGetUsageName(IOHIDElementGetUsagePage(currentElement), IOHIDElementGetUsage(currentElement), usageName);
        PsychSetStructArrayStringElement("usageName",		elementIndex, 	usageName,	 		elementStruct);

        PsychSetStructArrayDoubleElement("dataSize",		elementIndex, 	(double)IOHIDElementGetReportSize(currentElement), 	elementStruct);
        PsychSetStructArrayDoubleElement("rangeMin",		elementIndex, 	(double)IOHIDElementGetLogicalMin(currentElement), 	elementStruct);
        PsychSetStructArrayDoubleElement("rangeMax",		elementIndex, 	(double)IOHIDElementGetLogicalMax(currentElement), 	elementStruct);
        PsychSetStructArrayDoubleElement("scaledRangeMin",	elementIndex, 	(double)IOHIDElementGetPhysicalMin(currentElement), 	elementStruct);
        PsychSetStructArrayDoubleElement("scaledRangeMax",	elementIndex, 	(double)IOHIDElementGetPhysicalMax(currentElement),	elementStruct);
        PsychSetStructArrayDoubleElement("relative",		elementIndex, 	(double)IOHIDElementIsRelative(currentElement),	elementStruct);	//psych_bool flag
        PsychSetStructArrayDoubleElement("wrapping",		elementIndex, 	(double)IOHIDElementIsWrapping(currentElement),	elementStruct);	//psych_bool flag	
        PsychSetStructArrayDoubleElement("nonLinear",		elementIndex, 	(double)IOHIDElementIsNonLinear(currentElement),	elementStruct);	//psych_bool flag
        PsychSetStructArrayDoubleElement("preferredState",	elementIndex, 	(double)IOHIDElementHasPreferredState(currentElement),	elementStruct);	//psych_bool flag
        PsychSetStructArrayDoubleElement("nullState",		elementIndex, 	(double)IOHIDElementHasNullState(currentElement),	elementStruct);	//psych_bool flag
        PsychSetStructArrayDoubleElement("calMin",          elementIndex, 	(double)IOHIDElement_GetCalibrationMin(currentElement),	elementStruct);
        PsychSetStructArrayDoubleElement("calMax",          elementIndex, 	(double)IOHIDElement_GetCalibrationMax(currentElement),	elementStruct);
        PsychSetStructArrayDoubleElement("scalingMin",		elementIndex, 	(double) 0,     elementStruct);
        PsychSetStructArrayDoubleElement("scalingMax",		elementIndex, 	(double) 255,	elementStruct);
#else
        typeMask=HIDConvertElementTypeToMask (currentElement->type);
        PsychHIDGetTypeMaskStringFromTypeMask(typeMask, &typeMaskName);
        PsychSetStructArrayStringElement("typeMaskName",	elementIndex, 	typeMaskName,	 			elementStruct);
        PsychSetStructArrayStringElement("name",            elementIndex, 	currentElement->name,	 		elementStruct);
        PsychSetStructArrayDoubleElement("deviceIndex",		elementIndex, 	(double)deviceIndex, 			elementStruct);
        PsychSetStructArrayDoubleElement("elementIndex",	elementIndex, 	(double)elementIndex+1, 		elementStruct);
        PsychSetStructArrayDoubleElement("typeValue",		elementIndex, 	(double)currentElement->type, 		elementStruct);
        HIDGetTypeName(currentElement->type, elementTypeName);
        PsychSetStructArrayStringElement("typeName",		elementIndex, 	elementTypeName,	 		elementStruct);
        PsychSetStructArrayDoubleElement("usagePageValue",	elementIndex, 	(double)currentElement->usagePage, 	elementStruct);
        PsychSetStructArrayDoubleElement("usageValue",		elementIndex, 	(double)currentElement->usage, 		elementStruct);
        HIDGetUsageName (currentElement->usagePage, currentElement->usage, usageName);
        PsychSetStructArrayStringElement("usageName",		elementIndex, 	usageName,	 			elementStruct);
        PsychSetStructArrayDoubleElement("dataSize",		elementIndex, 	(double)currentElement->size, 		elementStruct);
        PsychSetStructArrayDoubleElement("rangeMin",		elementIndex, 	(double)currentElement->min, 		elementStruct);
        PsychSetStructArrayDoubleElement("rangeMax",		elementIndex, 	(double)currentElement->max, 		elementStruct);
        PsychSetStructArrayDoubleElement("scaledRangeMin",	elementIndex, 	(double)currentElement->scaledMin, 	elementStruct);
        PsychSetStructArrayDoubleElement("scaledRangeMax",	elementIndex, 	(double)currentElement->scaledMax,	elementStruct);
        PsychSetStructArrayDoubleElement("relative",		elementIndex, 	(double)currentElement->relative,	elementStruct);	//psych_bool flag
        PsychSetStructArrayDoubleElement("wrapping",		elementIndex, 	(double)currentElement->wrapping,	elementStruct);	//psych_bool flag	
        PsychSetStructArrayDoubleElement("nonLinear",		elementIndex, 	(double)currentElement->nonLinear,	elementStruct);	//psych_bool flag
        PsychSetStructArrayDoubleElement("preferredState",	elementIndex, 	(double)currentElement->preferredState,	elementStruct);	//psych_bool flag
        PsychSetStructArrayDoubleElement("nullState",		elementIndex, 	(double)currentElement->nullState,	elementStruct);	//psych_bool flag

        PsychSetStructArrayDoubleElement("calMin",		elementIndex, 	(double)currentElement->calMin,		elementStruct);	
        PsychSetStructArrayDoubleElement("calMax",		elementIndex, 	(double)currentElement->calMax,		elementStruct);	
        PsychSetStructArrayDoubleElement("scalingMin",		elementIndex, 	(double)currentElement->userMin,	elementStruct);	
        PsychSetStructArrayDoubleElement("scalingMax",		elementIndex, 	(double)currentElement->userMax,	elementStruct);	
#endif
        ++elementIndex; 
    }

    return(PsychError_none);	
}

#endif
