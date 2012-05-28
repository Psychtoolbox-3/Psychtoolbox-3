/*
  PsychToolbox3/Source/Common/PsychHID/PsychHIDGetDevices.c		
  
  PROJECTS: PsychHID
  
  PLATFORMS:  All  
  
  AUTHORS:
  Allen.Ingling@nyu.edu             awi 
  mario.kleiner@tuebingen.mpg.de    mk
      
  HISTORY:
  4/29/03  awi		Created.
  
  TO DO:
  
*/

#include "PsychHID.h"

static char useString[]= "devices=PsychHID('Devices' [, deviceClass])";
static char synopsisString[] =  "Return a struct array describing each connected USB HID device.\n"
				"'deviceClass' optionally selects for the class of input device. "
				"This is not supported on all operating systems and will be silently "
				"ignored if unsupported. On Linux you can select the following classes "
				"of input devices: 1 = MasterPointer, 2 = MasterKeyboard, 3 = SlavePointer "
				"4 = SlaveKeyboard, 5 = Floating slave device.\n\n"
                "Not all device properties are returned on all operating systems. A zero, "
                "empty or -1 value for a property in the returned structs can mean that "
                "the information could not be returned.\n";
static char seeAlsoString[] = "";

PsychError PSYCHHIDGetDevices(void) 
{
    pRecDevice 			currentDevice=NULL;
    
    const char *deviceFieldNames[]={"usagePageValue", "usageValue", "usageName", "index", "transport", "vendorID", "productID", "version", 
                                    "manufacturer", "product", "serialNumber", "locationID", "interfaceID", "totalElements", "features", "inputs", 
                                    "outputs", "collections", "axes", "buttons", "hats", "sliders", "dials", "wheels"};
    int numDeviceStructElements, numDeviceStructFieldNames=24, deviceIndex, deviceClass;
    PsychGenericScriptType	*deviceStruct;		
    char usageName[PSYCH_HID_MAX_DEVICE_ELEMENT_USAGE_NAME_LENGTH];
    char tmpString[1024];

    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));

    if (PsychCopyInIntegerArg(1, FALSE, &deviceClass)) {
        // Operating system specific enumeration of devices, selected by deviceClass:
        // Currently unsupported on OSX:
        #if PSYCH_SYSTEM != PSYCH_OSX
        return(PsychHIDEnumerateHIDInputDevices(deviceClass));
        #endif
    }

#if (PSYCH_SYSTEM == PSYCH_OSX) && defined(__LP64__)
    // Try to load all bundles from Psychtoolbox/PsychHardware/
    // This loads the HID_Utilities.framework bundle if it is present. The whole point of it is
    // to allow our statically compiled-in version of the library to find the location of
    // the XML file with the database of (vendorId, productId) -> (VendorName, ProductName) and
    // (usagePage, usage) -> (usageName) mappings.
    //
    // In practice, the XML file only serves as a fallback, and one that doesn't contain much
    // useful info for mainstream products, only for a few niche products. Given its limited
    // value, i think we can refrain from shipping the framework as part of Psychtoolbox and
    // just provide the option to use it (== its XML file) if users decide to install it themselves.
    
    sprintf(tmpString, "%sPsychHardware/", PsychRuntimeGetPsychtoolboxRoot(FALSE));
    CFStringRef urlString = CFStringCreateWithCString(kCFAllocatorDefault, tmpString, kCFStringEncodingASCII);
    CFURLRef directoryURL = CFURLCreateWithString(kCFAllocatorDefault, urlString, NULL);
    CFRelease(urlString);
    CFArrayRef bundleArray = CFBundleCreateBundlesFromDirectory(kCFAllocatorDefault, directoryURL, NULL);
    CFRelease(directoryURL);
    CFRelease(bundleArray);
    
    // This would detect if the bundle could be loaded:
    // CFBundleRef tCFBundleRef = CFBundleGetBundleWithIdentifier(CFSTR("com.apple.HID_Utilities"));
    // if (tCFBundleRef) printf("YEEEEEEEAAAAAAHHHH!\n");
#endif

    PsychHIDVerifyInit();
    numDeviceStructElements=(int)HIDCountDevices();
    PsychAllocOutStructArray(1, FALSE, numDeviceStructElements, numDeviceStructFieldNames, deviceFieldNames, &deviceStruct);
    deviceIndex=0;
    for(currentDevice=HIDGetFirstDevice(); currentDevice != NULL; currentDevice=HIDGetNextDevice(currentDevice)){

        // Code path for Linux, Windows and 32-Bit OSX 10.4 and later:
        #if (PSYCH_SYSTEM != PSYCH_OSX) || !defined(__LP64__)
            PsychSetStructArrayDoubleElement("usagePageValue",	deviceIndex, 	(double)currentDevice->usagePage,	deviceStruct);
            PsychSetStructArrayDoubleElement("usageValue",		deviceIndex, 	(double)currentDevice->usage, 		deviceStruct);

            PsychSetStructArrayStringElement("transport",		deviceIndex, 	currentDevice->transport, 		deviceStruct);
            PsychSetStructArrayDoubleElement("vendorID",		deviceIndex, 	(double)currentDevice->vendorID, 	deviceStruct);
            PsychSetStructArrayDoubleElement("productID",		deviceIndex, 	(double)currentDevice->productID, 	deviceStruct);
            PsychSetStructArrayDoubleElement("version",         deviceIndex, 	(double)currentDevice->version, 	deviceStruct);
            PsychSetStructArrayStringElement("manufacturer",	deviceIndex, 	currentDevice->manufacturer, 		deviceStruct);
            PsychSetStructArrayStringElement("product",         deviceIndex, 	currentDevice->product, 		deviceStruct);
            PsychSetStructArrayStringElement("serialNumber",	deviceIndex, 	currentDevice->serial, 			deviceStruct);
            PsychSetStructArrayDoubleElement("locationID",		deviceIndex, 	(double)currentDevice->locID, 		deviceStruct);

            PsychSetStructArrayDoubleElement("totalElements",	deviceIndex, 	(double)currentDevice->totalElements, 	deviceStruct);
            PsychSetStructArrayDoubleElement("features",		deviceIndex, 	(double)currentDevice->features, 	deviceStruct);
            PsychSetStructArrayDoubleElement("inputs",          deviceIndex, 	(double)currentDevice->inputs, 		deviceStruct);
            PsychSetStructArrayDoubleElement("outputs",         deviceIndex, 	(double)currentDevice->outputs, 	deviceStruct);
            PsychSetStructArrayDoubleElement("collections",		deviceIndex, 	(double)currentDevice->collections, 	deviceStruct);
            PsychSetStructArrayDoubleElement("axes",            deviceIndex, 	(double)currentDevice->axis, 		deviceStruct);
            PsychSetStructArrayDoubleElement("buttons",         deviceIndex, 	(double)currentDevice->buttons, 	deviceStruct);
            PsychSetStructArrayDoubleElement("hats",            deviceIndex, 	(double)currentDevice->hats, 		deviceStruct);
            PsychSetStructArrayDoubleElement("sliders",         deviceIndex, 	(double)currentDevice->sliders, 	deviceStruct);
            PsychSetStructArrayDoubleElement("dials",           deviceIndex, 	(double)currentDevice->dials, 		deviceStruct);
            PsychSetStructArrayDoubleElement("wheels",          deviceIndex, 	(double)currentDevice->wheels, 		deviceStruct);
        #endif

        // OSX specific:
        #if PSYCH_SYSTEM == PSYCH_OSX
            // OSX 64-Bit specific: OSX 10.5 and later:
            #ifdef __LP64__
                PsychSetStructArrayDoubleElement("usagePageValue", deviceIndex, (double) IOHIDDevice_GetPrimaryUsagePage(currentDevice), deviceStruct);
                PsychSetStructArrayDoubleElement("usageValue", deviceIndex, (double) IOHIDDevice_GetPrimaryUsage(currentDevice), deviceStruct);

                sprintf(usageName, "");
                CFStringRef cfusageName = HIDCopyUsageName(IOHIDDevice_GetPrimaryUsagePage(currentDevice), IOHIDDevice_GetPrimaryUsage(currentDevice));
                if (cfusageName && (CFStringGetLength(cfusageName) > 0)) {
                    CFStringGetCString(cfusageName, usageName, sizeof(usageName), kCFStringEncodingASCII);
                    CFRelease(cfusageName);
                } else {
                    // Use our fallback from HID Utilties v1.0, implemented below this function at the end of the file:
                    HIDGetUsageName(IOHIDDevice_GetPrimaryUsagePage(currentDevice), IOHIDDevice_GetPrimaryUsage(currentDevice), usageName);
                }

                PsychSetStructArrayDoubleElement("vendorID", deviceIndex, (double) IOHIDDevice_GetVendorID(currentDevice), deviceStruct);
                PsychSetStructArrayDoubleElement("productID", deviceIndex, (double)IOHIDDevice_GetProductID(currentDevice), deviceStruct);
                PsychSetStructArrayDoubleElement("version", deviceIndex, (double) IOHIDDevice_GetVersionNumber(currentDevice), deviceStruct);
                PsychSetStructArrayDoubleElement("locationID", deviceIndex, (double) IOHIDDevice_GetLocationID(currentDevice), deviceStruct);

                sprintf(tmpString, "");
                cfusageName = IOHIDDevice_GetManufacturer(currentDevice);
                if (cfusageName && (CFStringGetLength(cfusageName) > 0)) {
                    CFStringGetCString(cfusageName, tmpString, sizeof(tmpString), kCFStringEncodingASCII);
                } else {
                    HIDGetVendorNameFromVendorID(IOHIDDevice_GetVendorID(currentDevice), tmpString);
                }
                PsychSetStructArrayStringElement("manufacturer", deviceIndex, tmpString, deviceStruct);

                sprintf(tmpString, "");
                cfusageName = IOHIDDevice_GetProduct(currentDevice);
                if (cfusageName && (CFStringGetLength(cfusageName) > 0)) {
                    CFStringGetCString(cfusageName, tmpString, sizeof(tmpString), kCFStringEncodingASCII);
                } else {
                    HIDGetProductNameFromVendorProductID(IOHIDDevice_GetVendorID(currentDevice), IOHIDDevice_GetProductID(currentDevice), tmpString);
                }
                PsychSetStructArrayStringElement("product", deviceIndex, tmpString, deviceStruct);
        
                sprintf(tmpString, "");
                cfusageName = IOHIDDevice_GetTransport(currentDevice);
                if (cfusageName && (CFStringGetLength(cfusageName) > 0)) {
                    CFStringGetCString(cfusageName, tmpString, sizeof(tmpString), kCFStringEncodingASCII);
                }
                PsychSetStructArrayStringElement("transport", deviceIndex, tmpString, deviceStruct);

                sprintf(tmpString, "");
                cfusageName = IOHIDDevice_GetSerialNumber(currentDevice);
                if (cfusageName && (CFStringGetLength(cfusageName) > 0)) {
                    CFStringGetCString(cfusageName, tmpString, sizeof(tmpString), kCFStringEncodingASCII);
                }
                PsychSetStructArrayStringElement("serialNumber", deviceIndex, tmpString, deviceStruct);

                PsychSetStructArrayDoubleElement("collections",     deviceIndex,    HIDCountDeviceElements(currentDevice, kHIDElementTypeCollection), deviceStruct);
                PsychSetStructArrayDoubleElement("totalElements",   deviceIndex,    HIDCountDeviceElements(currentDevice, kHIDElementTypeIO), deviceStruct);
                PsychSetStructArrayDoubleElement("features",   deviceIndex,    HIDCountDeviceElements(currentDevice, kHIDElementTypeFeature), deviceStruct);
                PsychSetStructArrayDoubleElement("inputs",   deviceIndex,    HIDCountDeviceElements(currentDevice, kHIDElementTypeInput), deviceStruct);
                PsychSetStructArrayDoubleElement("outputs",   deviceIndex,    HIDCountDeviceElements(currentDevice, kHIDElementTypeOutput), deviceStruct);
    
                // Iterate over all device input elements and count buttons, sliders, axis, hats etc.:
                pRecElement currentElement, lastElement = NULL;
                long usagePage, usage;
                unsigned int axis = 0, sliders = 0, dials = 0, wheels = 0, hats = 0, buttons = 0;
                for(currentElement = HIDGetFirstDeviceElement(currentDevice, kHIDElementTypeInput); 
                    (currentElement != NULL) && (currentElement != lastElement);
                    currentElement = HIDGetNextDeviceElement(currentElement, kHIDElementTypeInput)) {
                    lastElement = currentElement;
                    
                    usagePage = IOHIDElementGetUsagePage(currentElement);
                    usage = IOHIDElementGetUsage(currentElement);
                    
                    switch (usagePage) {
                        case kHIDPage_GenericDesktop:
                            switch (usage) {
                                case kHIDUsage_GD_X:
                                case kHIDUsage_GD_Y:
                                case kHIDUsage_GD_Z:
                                case kHIDUsage_GD_Rx:
                                case kHIDUsage_GD_Ry:
                                case kHIDUsage_GD_Rz:
                                    axis++;
                                break;
                                    
                                case kHIDUsage_GD_Slider:
                                    sliders++;
                                break;

                                case kHIDUsage_GD_Dial:
                                    dials++;
                                    break;

                                case kHIDUsage_GD_Wheel:
                                    wheels++;
                                    break;

                                case kHIDUsage_GD_Hatswitch:
                                    hats++;
                                    break;
                                default:
                                    break;
                            }
                        break;
                            
                        case kHIDPage_Button:
                            buttons++;
                        break;
                            
                        default:
                        break;
                    }
                }

                PsychSetStructArrayDoubleElement("axes",            deviceIndex, 	(double) axis, 		deviceStruct);
                PsychSetStructArrayDoubleElement("buttons",         deviceIndex, 	(double) buttons, 	deviceStruct);
                PsychSetStructArrayDoubleElement("hats",            deviceIndex, 	(double) hats, 		deviceStruct);
                PsychSetStructArrayDoubleElement("sliders",         deviceIndex, 	(double) sliders, 	deviceStruct);
                PsychSetStructArrayDoubleElement("dials",           deviceIndex, 	(double) dials, 		deviceStruct);
                PsychSetStructArrayDoubleElement("wheels",          deviceIndex, 	(double) wheels, 		deviceStruct);
            #else
                // 32-Bit OSX 10.4 and later:
                HIDGetUsageName(currentDevice->usagePage, currentDevice->usage, usageName);
            #endif
        
            // Store dummy value -1 to mark interfaceID as invalid/unknown on OSX:
            PsychSetStructArrayDoubleElement("interfaceID",	deviceIndex, (double) -1, deviceStruct);
        #else
            // Linux, Windows:

            // USB interface id only available on non OS/X:
            PsychSetStructArrayDoubleElement("interfaceID",	deviceIndex, (double)currentDevice->interfaceId, deviceStruct);

            // TODO FIXME Usage name: Mapping of usagePage + usage to human readable string
            // is to be done for Linux/Windows: HIDGetUsageName (currentDevice->usagePage, currentDevice->usage, usageName);
            sprintf(usageName, "TBD");        
        #endif

        // UsageName as parsed in os specific code above:
        PsychSetStructArrayStringElement("usageName", deviceIndex, usageName, deviceStruct);
        
        // Logical device index for ptb:
        PsychSetStructArrayDoubleElement("index", deviceIndex, (double) deviceIndex+1, deviceStruct);

        // Next device...
        ++deviceIndex; 
    }

    return(PsychError_none);	
}

#if (PSYCH_SYSTEM == PSYCH_OSX) && defined(__LP64__)
// This routine is a transplanted and modified version from Apple HID Utiltities v1.0, as allowed by
// its license, which would be BSD-3 style for unmodified distribution of the full package, but has
// no restrictions or requirements at all if only parts are copied, or code is modified after copy,
// as here, where we only extract one routine and then modify it.
//
// For reference:
// The source code of the original HID Utilties can be found in PsychSourceGL/Cohorts/HID Utilities/
// this is from file HID_Utilities.c
//
// ---------------------------------
// returns C string usage given usage page and usage passed in as parameters (see IOUSBHIDParser.h)
// returns usage page and usage values in string form for unknown values

void HIDGetUsageName (const long valueUsagePage, const long valueUsage, char * cstrName)
{
    switch (valueUsagePage)
    {
        case kHIDPage_Undefined:
			switch (valueUsage)
        {
            default: sprintf (cstrName, "Undefined Page, Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_GenericDesktop:
            switch (valueUsage)
        {
            case kHIDUsage_GD_Pointer: sprintf (cstrName, "Pointer"); break;
            case kHIDUsage_GD_Mouse: sprintf (cstrName, "Mouse"); break;
            case kHIDUsage_GD_Joystick: sprintf (cstrName, "Joystick"); break;
            case kHIDUsage_GD_GamePad: sprintf (cstrName, "GamePad"); break;
            case kHIDUsage_GD_Keyboard: sprintf (cstrName, "Keyboard"); break;
            case kHIDUsage_GD_Keypad: sprintf (cstrName, "Keypad"); break;
            case kHIDUsage_GD_MultiAxisController: sprintf (cstrName, "Multi-Axis Controller"); break;
                
            case kHIDUsage_GD_X: sprintf (cstrName, "X-Axis"); break;
            case kHIDUsage_GD_Y: sprintf (cstrName, "Y-Axis"); break;
            case kHIDUsage_GD_Z: sprintf (cstrName, "Z-Axis"); break;
            case kHIDUsage_GD_Rx: sprintf (cstrName, "X-Rotation"); break;
            case kHIDUsage_GD_Ry: sprintf (cstrName, "Y-Rotation"); break;
            case kHIDUsage_GD_Rz: sprintf (cstrName, "Z-Rotation"); break;
            case kHIDUsage_GD_Slider: sprintf (cstrName, "Slider"); break;
            case kHIDUsage_GD_Dial: sprintf (cstrName, "Dial"); break;
            case kHIDUsage_GD_Wheel: sprintf (cstrName, "Wheel"); break;
            case kHIDUsage_GD_Hatswitch: sprintf (cstrName, "Hatswitch"); break;
            case kHIDUsage_GD_CountedBuffer: sprintf (cstrName, "Counted Buffer"); break;
            case kHIDUsage_GD_ByteCount: sprintf (cstrName, "Byte Count"); break;
            case kHIDUsage_GD_MotionWakeup: sprintf (cstrName, "Motion Wakeup"); break;
            case kHIDUsage_GD_Start: sprintf (cstrName, "Start"); break;
            case kHIDUsage_GD_Select: sprintf (cstrName, "Select"); break;
                
            case kHIDUsage_GD_Vx: sprintf (cstrName, "X-Velocity"); break;
            case kHIDUsage_GD_Vy: sprintf (cstrName, "Y-Velocity"); break;
            case kHIDUsage_GD_Vz: sprintf (cstrName, "Z-Velocity"); break;
            case kHIDUsage_GD_Vbrx: sprintf (cstrName, "X-Rotation Velocity"); break;
            case kHIDUsage_GD_Vbry: sprintf (cstrName, "Y-Rotation Velocity"); break;
            case kHIDUsage_GD_Vbrz: sprintf (cstrName, "Z-Rotation Velocity"); break;
            case kHIDUsage_GD_Vno: sprintf (cstrName, "Vno"); break;
                
            case kHIDUsage_GD_SystemControl: sprintf (cstrName, "System Control"); break;
            case kHIDUsage_GD_SystemPowerDown: sprintf (cstrName, "System Power Down"); break;
            case kHIDUsage_GD_SystemSleep: sprintf (cstrName, "System Sleep"); break;
            case kHIDUsage_GD_SystemWakeUp: sprintf (cstrName, "System Wake Up"); break;
            case kHIDUsage_GD_SystemContextMenu: sprintf (cstrName, "System Context Menu"); break;
            case kHIDUsage_GD_SystemMainMenu: sprintf (cstrName, "System Main Menu"); break;
            case kHIDUsage_GD_SystemAppMenu: sprintf (cstrName, "System App Menu"); break;
            case kHIDUsage_GD_SystemMenuHelp: sprintf (cstrName, "System Menu Help"); break;
            case kHIDUsage_GD_SystemMenuExit: sprintf (cstrName, "System Menu Exit"); break;
            case kHIDUsage_GD_SystemMenu: sprintf (cstrName, "System Menu"); break;
            case kHIDUsage_GD_SystemMenuRight: sprintf (cstrName, "System Menu Right"); break;
            case kHIDUsage_GD_SystemMenuLeft: sprintf (cstrName, "System Menu Left"); break;
            case kHIDUsage_GD_SystemMenuUp: sprintf (cstrName, "System Menu Up"); break;
            case kHIDUsage_GD_SystemMenuDown: sprintf (cstrName, "System Menu Down"); break;
                
            case kHIDUsage_GD_DPadUp: sprintf (cstrName, "DPad Up"); break;
            case kHIDUsage_GD_DPadDown: sprintf (cstrName, "DPad Down"); break;
            case kHIDUsage_GD_DPadRight: sprintf (cstrName, "DPad Right"); break;
            case kHIDUsage_GD_DPadLeft: sprintf (cstrName, "DPad Left"); break;
                
            case kHIDUsage_GD_Reserved: sprintf (cstrName, "Reserved"); break;
                
            default: sprintf (cstrName, "Generic Desktop Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Simulation:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Simulation Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_VR:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "VR Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Sport:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Sport Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Game:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Game Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_KeyboardOrKeypad:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Keyboard Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_LEDs:
            switch (valueUsage)
        {
				// some LED usages
            case kHIDUsage_LED_IndicatorRed: sprintf (cstrName, "Red LED"); break;
            case kHIDUsage_LED_IndicatorGreen: sprintf (cstrName, "Green LED"); break;
            case kHIDUsage_LED_IndicatorAmber: sprintf (cstrName, "Amber LED"); break;
            case kHIDUsage_LED_GenericIndicator: sprintf (cstrName, "Generic LED"); break;
            case kHIDUsage_LED_SystemSuspend: sprintf (cstrName, "System Suspend LED"); break;
            case kHIDUsage_LED_ExternalPowerConnected: sprintf (cstrName, "External Power LED"); break;
            default: sprintf (cstrName, "LED Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Button:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Button #%ld", valueUsage); break;
        }
            break;
        case kHIDPage_Ordinal:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Ordinal Instance %lx", valueUsage); break;
        }
            break;
        case kHIDPage_Telephony:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Telephony Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Consumer:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Consumer Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Digitizer:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Digitizer Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_PID:
			if (((valueUsage >= 0x02) && (valueUsage <= 0x1F)) || ((valueUsage >= 0x29) && (valueUsage <= 0x2F)) ||
                ((valueUsage >= 0x35) && (valueUsage <= 0x3F)) || ((valueUsage >= 0x44) && (valueUsage <= 0x4F)) ||
                (valueUsage == 0x8A) || (valueUsage == 0x93)  || ((valueUsage >= 0x9D) && (valueUsage <= 0x9E)) ||
                ((valueUsage >= 0xA1) && (valueUsage <= 0xA3)) || ((valueUsage >= 0xAD) && (valueUsage <= 0xFFFF)))
                sprintf (cstrName, "PID Reserved");
			else
				switch (valueUsage)
            {
                case 0x00: sprintf (cstrName, "PID Undefined Usage"); break;
                case kHIDUsage_PID_PhysicalInterfaceDevice: sprintf (cstrName, "Physical Interface Device"); break;
                case kHIDUsage_PID_Normal: sprintf (cstrName, "Normal Force"); break;
                    
                case kHIDUsage_PID_SetEffectReport: sprintf (cstrName, "Set Effect Report"); break;
                case kHIDUsage_PID_EffectBlockIndex: sprintf (cstrName, "Effect Block Index"); break;
                case kHIDUsage_PID_ParamBlockOffset: sprintf (cstrName, "Parameter Block Offset"); break;
                case kHIDUsage_PID_ROM_Flag: sprintf (cstrName, "ROM Flag"); break;
                    
                case kHIDUsage_PID_EffectType: sprintf (cstrName, "Effect Type"); break;
                case kHIDUsage_PID_ET_ConstantForce: sprintf (cstrName, "Effect Type Constant Force"); break;
                case kHIDUsage_PID_ET_Ramp: sprintf (cstrName, "Effect Type Ramp"); break;
                case kHIDUsage_PID_ET_CustomForceData: sprintf (cstrName, "Effect Type Custom Force Data"); break;
                case kHIDUsage_PID_ET_Square: sprintf (cstrName, "Effect Type Square"); break;
                case kHIDUsage_PID_ET_Sine: sprintf (cstrName, "Effect Type Sine"); break;
                case kHIDUsage_PID_ET_Triangle: sprintf (cstrName, "Effect Type Triangle"); break;
                case kHIDUsage_PID_ET_SawtoothUp: sprintf (cstrName, "Effect Type Sawtooth Up"); break;
                case kHIDUsage_PID_ET_SawtoothDown: sprintf (cstrName, "Effect Type Sawtooth Down"); break;
                case kHIDUsage_PID_ET_Spring: sprintf (cstrName, "Effect Type Spring"); break;
                case kHIDUsage_PID_ET_Damper: sprintf (cstrName, "Effect Type Damper"); break;
                case kHIDUsage_PID_ET_Inertia: sprintf (cstrName, "Effect Type Inertia"); break;
                case kHIDUsage_PID_ET_Friction: sprintf (cstrName, "Effect Type Friction"); break;
                case kHIDUsage_PID_Duration: sprintf (cstrName, "Effect Duration"); break;
                case kHIDUsage_PID_SamplePeriod: sprintf (cstrName, "Effect Sample Period"); break;
                case kHIDUsage_PID_Gain: sprintf (cstrName, "Effect Gain"); break;
                case kHIDUsage_PID_TriggerButton: sprintf (cstrName, "Effect Trigger Button"); break;
                case kHIDUsage_PID_TriggerRepeatInterval: sprintf (cstrName, "Effect Trigger Repeat Interval"); break;
                    
                case kHIDUsage_PID_AxesEnable: sprintf (cstrName, "Axis Enable"); break;
                case kHIDUsage_PID_DirectionEnable: sprintf (cstrName, "Direction Enable"); break;
                    
                case kHIDUsage_PID_Direction: sprintf (cstrName, "Direction"); break;
                    
                case kHIDUsage_PID_TypeSpecificBlockOffset: sprintf (cstrName, "Type Specific Block Offset"); break;
                    
                case kHIDUsage_PID_BlockType: sprintf (cstrName, "Block Type"); break;
                    
                case kHIDUsage_PID_SetEnvelopeReport: sprintf (cstrName, "Set Envelope Report"); break;
                case kHIDUsage_PID_AttackLevel: sprintf (cstrName, "Envelope Attack Level"); break;
                case kHIDUsage_PID_AttackTime: sprintf (cstrName, "Envelope Attack Time"); break;
                case kHIDUsage_PID_FadeLevel: sprintf (cstrName, "Envelope Fade Level"); break;
                case kHIDUsage_PID_FadeTime: sprintf (cstrName, "Envelope Fade Time"); break;
                    
                case kHIDUsage_PID_SetConditionReport: sprintf (cstrName, "Set Condition Report"); break;
                case kHIDUsage_PID_CP_Offset: sprintf (cstrName, "Condition CP Offset"); break;
                case kHIDUsage_PID_PositiveCoefficient: sprintf (cstrName, "Condition Positive Coefficient"); break;
                case kHIDUsage_PID_NegativeCoefficient: sprintf (cstrName, "Condition Negative Coefficient"); break;
                case kHIDUsage_PID_PositiveSaturation: sprintf (cstrName, "Condition Positive Saturation"); break;
                case kHIDUsage_PID_NegativeSaturation: sprintf (cstrName, "Condition Negative Saturation"); break;
                case kHIDUsage_PID_DeadBand: sprintf (cstrName, "Condition Dead Band"); break;
                    
                case kHIDUsage_PID_DownloadForceSample: sprintf (cstrName, "Download Force Sample"); break;
                case kHIDUsage_PID_IsochCustomForceEnable: sprintf (cstrName, "Isoch Custom Force Enable"); break;
                    
                case kHIDUsage_PID_CustomForceDataReport: sprintf (cstrName, "Custom Force Data Report"); break;
                case kHIDUsage_PID_CustomForceData: sprintf (cstrName, "Custom Force Data"); break;
                    
                case kHIDUsage_PID_CustomForceVendorDefinedData: sprintf (cstrName, "Custom Force Vendor Defined Data"); break;
                case kHIDUsage_PID_SetCustomForceReport: sprintf (cstrName, "Set Custom Force Report"); break;
                case kHIDUsage_PID_CustomForceDataOffset: sprintf (cstrName, "Custom Force Data Offset"); break;
                case kHIDUsage_PID_SampleCount: sprintf (cstrName, "Custom Force Sample Count"); break;
                    
                case kHIDUsage_PID_SetPeriodicReport: sprintf (cstrName, "Set Periodic Report"); break;
                case kHIDUsage_PID_Offset: sprintf (cstrName, "Periodic Offset"); break;
                case kHIDUsage_PID_Magnitude: sprintf (cstrName, "Periodic Magnitude"); break;
                case kHIDUsage_PID_Phase: sprintf (cstrName, "Periodic Phase"); break;
                case kHIDUsage_PID_Period: sprintf (cstrName, "Periodic Period"); break;
                    
                case kHIDUsage_PID_SetConstantForceReport: sprintf (cstrName, "Set Constant Force Report"); break;
                    
                case kHIDUsage_PID_SetRampForceReport: sprintf (cstrName, "Set Ramp Force Report"); break;
                case kHIDUsage_PID_RampStart: sprintf (cstrName, "Ramp Start"); break;
                case kHIDUsage_PID_RampEnd: sprintf (cstrName, "Ramp End"); break;
                    
                case kHIDUsage_PID_EffectOperationReport: sprintf (cstrName, "Effect Operation Report"); break;
                    
                case kHIDUsage_PID_EffectOperation: sprintf (cstrName, "Effect Operation"); break;
                case kHIDUsage_PID_OpEffectStart: sprintf (cstrName, "Op Effect Start"); break;
                case kHIDUsage_PID_OpEffectStartSolo: sprintf (cstrName, "Op Effect Start Solo"); break;
                case kHIDUsage_PID_OpEffectStop: sprintf (cstrName, "Op Effect Stop"); break;
                case kHIDUsage_PID_LoopCount: sprintf (cstrName, "Op Effect Loop Count"); break;
                    
                case kHIDUsage_PID_DeviceGainReport: sprintf (cstrName, "Device Gain Report"); break;
                case kHIDUsage_PID_DeviceGain: sprintf (cstrName, "Device Gain"); break;
                    
                case kHIDUsage_PID_PoolReport: sprintf (cstrName, "PID Pool Report"); break;
                case kHIDUsage_PID_RAM_PoolSize: sprintf (cstrName, "RAM Pool Size"); break;
                case kHIDUsage_PID_ROM_PoolSize: sprintf (cstrName, "ROM Pool Size"); break;
                case kHIDUsage_PID_ROM_EffectBlockCount: sprintf (cstrName, "ROM Effect Block Count"); break;
                case kHIDUsage_PID_SimultaneousEffectsMax: sprintf (cstrName, "Simultaneous Effects Max"); break;
                case kHIDUsage_PID_PoolAlignment: sprintf (cstrName, "Pool Alignment"); break;
                    
                case kHIDUsage_PID_PoolMoveReport: sprintf (cstrName, "PID Pool Move Report"); break;
                case kHIDUsage_PID_MoveSource: sprintf (cstrName, "Move Source"); break;
                case kHIDUsage_PID_MoveDestination: sprintf (cstrName, "Move Destination"); break;
                case kHIDUsage_PID_MoveLength: sprintf (cstrName, "Move Length"); break;
                    
                case kHIDUsage_PID_BlockLoadReport: sprintf (cstrName, "PID Block Load Report"); break;
                    
                case kHIDUsage_PID_BlockLoadStatus: sprintf (cstrName, "Block Load Status"); break;
                case kHIDUsage_PID_BlockLoadSuccess: sprintf (cstrName, "Block Load Success"); break;
                case kHIDUsage_PID_BlockLoadFull: sprintf (cstrName, "Block Load Full"); break;
                case kHIDUsage_PID_BlockLoadError: sprintf (cstrName, "Block Load Error"); break;
                case kHIDUsage_PID_BlockHandle: sprintf (cstrName, "Block Handle"); break;
                    
                case kHIDUsage_PID_BlockFreeReport: sprintf (cstrName, "PID Block Free Report"); break;
                    
                case kHIDUsage_PID_TypeSpecificBlockHandle: sprintf (cstrName, "Type Specific Block Handle"); break;
                    
                case kHIDUsage_PID_StateReport: sprintf (cstrName, "PID State Report"); break;
                case kHIDUsage_PID_EffectPlaying: sprintf (cstrName, "Effect Playing"); break;
                    
                case kHIDUsage_PID_DeviceControlReport: sprintf (cstrName, "PID Device Control Report"); break;
                    
                case kHIDUsage_PID_DeviceControl: sprintf (cstrName, "PID Device Control"); break;
                case kHIDUsage_PID_DC_EnableActuators: sprintf (cstrName, "Device Control Enable Actuators"); break;
                case kHIDUsage_PID_DC_DisableActuators: sprintf (cstrName, "Device Control Disable Actuators"); break;
                case kHIDUsage_PID_DC_StopAllEffects: sprintf (cstrName, "Device Control Stop All Effects"); break;
                case kHIDUsage_PID_DC_DeviceReset: sprintf (cstrName, "Device Control Reset"); break;
                case kHIDUsage_PID_DC_DevicePause: sprintf (cstrName, "Device Control Pause"); break;
                case kHIDUsage_PID_DC_DeviceContinue: sprintf (cstrName, "Device Control Continue"); break;
                case kHIDUsage_PID_DevicePaused: sprintf (cstrName, "Device Paused"); break;
                case kHIDUsage_PID_ActuatorsEnabled: sprintf (cstrName, "Actuators Enabled"); break;
                case kHIDUsage_PID_SafetySwitch: sprintf (cstrName, "Safety Switch"); break;
                case kHIDUsage_PID_ActuatorOverrideSwitch: sprintf (cstrName, "Actuator Override Switch"); break;
                case kHIDUsage_PID_ActuatorPower: sprintf (cstrName, "Actuator Power"); break;
                case kHIDUsage_PID_StartDelay: sprintf (cstrName, "Start Delay"); break;
                    
                case kHIDUsage_PID_ParameterBlockSize: sprintf (cstrName, "Parameter Block Size"); break;
                case kHIDUsage_PID_DeviceManagedPool: sprintf (cstrName, "Device Managed Pool"); break;
                case kHIDUsage_PID_SharedParameterBlocks: sprintf (cstrName, "Shared Parameter Blocks"); break;
                    
                case kHIDUsage_PID_CreateNewEffectReport: sprintf (cstrName, "Create New Effect Report"); break;
                case kHIDUsage_PID_RAM_PoolAvailable: sprintf (cstrName, "RAM Pool Available"); break;
                default: sprintf (cstrName, "PID Usage 0x%lx", valueUsage); break;
            }
            break;
        case kHIDPage_Unicode:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Unicode Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_PowerDevice:
			if (((valueUsage >= 0x06) && (valueUsage <= 0x0F)) || ((valueUsage >= 0x26) && (valueUsage <= 0x2F)) ||
                ((valueUsage >= 0x39) && (valueUsage <= 0x3F)) || ((valueUsage >= 0x48) && (valueUsage <= 0x4F)) ||
                ((valueUsage >= 0x58) && (valueUsage <= 0x5F)) || (valueUsage == 0x6A) ||
                ((valueUsage >= 0x74) && (valueUsage <= 0xFC)))
                sprintf (cstrName, "Power Device Reserved");
			else
				switch (valueUsage)
            {
                case kHIDUsage_PD_Undefined: sprintf (cstrName, "Power Device Undefined Usage"); break;
                case kHIDUsage_PD_iName: sprintf (cstrName, "Power Device Name Index"); break;
                case kHIDUsage_PD_PresentStatus: sprintf (cstrName, "Power Device Present Status"); break;
                case kHIDUsage_PD_ChangedStatus: sprintf (cstrName, "Power Device Changed Status"); break;
                case kHIDUsage_PD_UPS: sprintf (cstrName, "Uninterruptible Power Supply"); break;
                case kHIDUsage_PD_PowerSupply: sprintf (cstrName, "Power Supply"); break;
                    
                case kHIDUsage_PD_BatterySystem: sprintf (cstrName, "Battery System Power Module"); break;
                case kHIDUsage_PD_BatterySystemID: sprintf (cstrName, "Battery System ID"); break;
                case kHIDUsage_PD_Battery: sprintf (cstrName, "Battery"); break;
                case kHIDUsage_PD_BatteryID: sprintf (cstrName, "Battery ID"); break;
                case kHIDUsage_PD_Charger: sprintf (cstrName, "Charger"); break;
                case kHIDUsage_PD_ChargerID: sprintf (cstrName, "Charger ID"); break;
                case kHIDUsage_PD_PowerConverter: sprintf (cstrName, "Power Converter Power Module"); break;
                case kHIDUsage_PD_PowerConverterID: sprintf (cstrName, "Power Converter ID"); break;
                case kHIDUsage_PD_OutletSystem: sprintf (cstrName, "Outlet System power module"); break;
                case kHIDUsage_PD_OutletSystemID: sprintf (cstrName, "Outlet System ID"); break;
                case kHIDUsage_PD_Input: sprintf (cstrName, "Power Device Input"); break;
                case kHIDUsage_PD_InputID: sprintf (cstrName, "Power Device Input ID"); break;
                case kHIDUsage_PD_Output: sprintf (cstrName, "Power Device Output"); break;
                case kHIDUsage_PD_OutputID: sprintf (cstrName, "Power Device Output ID"); break;
                case kHIDUsage_PD_Flow: sprintf (cstrName, "Power Device Flow"); break;
                case kHIDUsage_PD_FlowID: sprintf (cstrName, "Power Device Flow ID"); break;
                case kHIDUsage_PD_Outlet: sprintf (cstrName, "Power Device Outlet"); break;
                case kHIDUsage_PD_OutletID: sprintf (cstrName, "Power Device Outlet ID"); break;
                case kHIDUsage_PD_Gang: sprintf (cstrName, "Power Device Gang"); break;
                case kHIDUsage_PD_GangID: sprintf (cstrName, "Power Device Gang ID"); break;
                case kHIDUsage_PD_PowerSummary: sprintf (cstrName, "Power Device Power Summary"); break;
                case kHIDUsage_PD_PowerSummaryID: sprintf (cstrName, "Power Device Power Summary ID"); break;
                    
                case kHIDUsage_PD_Voltage: sprintf (cstrName, "Power Device Voltage"); break;
                case kHIDUsage_PD_Current: sprintf (cstrName, "Power Device Current"); break;
                case kHIDUsage_PD_Frequency: sprintf (cstrName, "Power Device Frequency"); break;
                case kHIDUsage_PD_ApparentPower: sprintf (cstrName, "Power Device Apparent Power"); break;
                case kHIDUsage_PD_ActivePower: sprintf (cstrName, "Power Device RMS Power"); break;
                case kHIDUsage_PD_PercentLoad: sprintf (cstrName, "Power Device Percent Load"); break;
                case kHIDUsage_PD_Temperature: sprintf (cstrName, "Power Device Temperature"); break;
                case kHIDUsage_PD_Humidity: sprintf (cstrName, "Power Device Humidity"); break;
                case kHIDUsage_PD_BadCount: sprintf (cstrName, "Power Device Bad Condition Count"); break;
                    
                case kHIDUsage_PD_ConfigVoltage: sprintf (cstrName, "Power Device Nominal Voltage"); break;
                case kHIDUsage_PD_ConfigCurrent: sprintf (cstrName, "Power Device Nominal Current"); break;
                case kHIDUsage_PD_ConfigFrequency: sprintf (cstrName, "Power Device Nominal Frequency"); break;
                case kHIDUsage_PD_ConfigApparentPower: sprintf (cstrName, "Power Device Nominal Apparent Power"); break;
                case kHIDUsage_PD_ConfigActivePower: sprintf (cstrName, "Power Device Nominal RMS Power"); break;
                case kHIDUsage_PD_ConfigPercentLoad: sprintf (cstrName, "Power Device Nominal Percent Load"); break;
                case kHIDUsage_PD_ConfigTemperature: sprintf (cstrName, "Power Device Nominal Temperature"); break;
                    
                case kHIDUsage_PD_ConfigHumidity: sprintf (cstrName, "Power Device Nominal Humidity"); break;
                case kHIDUsage_PD_SwitchOnControl: sprintf (cstrName, "Power Device Switch On Control"); break;
                case kHIDUsage_PD_SwitchOffControl: sprintf (cstrName, "Power Device Switch Off Control"); break;
                case kHIDUsage_PD_ToggleControl: sprintf (cstrName, "Power Device Toogle Sequence Control"); break;
                case kHIDUsage_PD_LowVoltageTransfer: sprintf (cstrName, "Power Device Min Transfer Voltage"); break;
                case kHIDUsage_PD_HighVoltageTransfer: sprintf (cstrName, "Power Device Max Transfer Voltage"); break;
                case kHIDUsage_PD_DelayBeforeReboot: sprintf (cstrName, "Power Device Delay Before Reboot"); break;
                case kHIDUsage_PD_DelayBeforeStartup: sprintf (cstrName, "Power Device Delay Before Startup"); break;
                case kHIDUsage_PD_DelayBeforeShutdown: sprintf (cstrName, "Power Device Delay Before Shutdown"); break;
                case kHIDUsage_PD_Test: sprintf (cstrName, "Power Device Test Request/Result"); break;
                case kHIDUsage_PD_ModuleReset: sprintf (cstrName, "Power Device Reset Request/Result"); break;
                case kHIDUsage_PD_AudibleAlarmControl: sprintf (cstrName, "Power Device Audible Alarm Control"); break;
                    
                case kHIDUsage_PD_Present: sprintf (cstrName, "Power Device Present"); break;
                case kHIDUsage_PD_Good: sprintf (cstrName, "Power Device Good"); break;
                case kHIDUsage_PD_InternalFailure: sprintf (cstrName, "Power Device Internal Failure"); break;
                case kHIDUsage_PD_VoltageOutOfRange: sprintf (cstrName, "Power Device Voltage Out Of Range"); break;
                case kHIDUsage_PD_FrequencyOutOfRange: sprintf (cstrName, "Power Device Frequency Out Of Range"); break;
                case kHIDUsage_PD_Overload: sprintf (cstrName, "Power Device Overload"); break;
                case kHIDUsage_PD_OverCharged: sprintf (cstrName, "Power Device Over Charged"); break;
                case kHIDUsage_PD_OverTemperature: sprintf (cstrName, "Power Device Over Temperature"); break;
                case kHIDUsage_PD_ShutdownRequested: sprintf (cstrName, "Power Device Shutdown Requested"); break;
                    
                case kHIDUsage_PD_ShutdownImminent: sprintf (cstrName, "Power Device Shutdown Imminent"); break;
                case kHIDUsage_PD_SwitchOnOff: sprintf (cstrName, "Power Device On/Off Switch Status"); break;
                case kHIDUsage_PD_Switchable: sprintf (cstrName, "Power Device Switchable"); break;
                case kHIDUsage_PD_Used: sprintf (cstrName, "Power Device Used"); break;
                case kHIDUsage_PD_Boost: sprintf (cstrName, "Power Device Boosted"); break;
                case kHIDUsage_PD_Buck: sprintf (cstrName, "Power Device Bucked"); break;
                case kHIDUsage_PD_Initialized: sprintf (cstrName, "Power Device Initialized"); break;
                case kHIDUsage_PD_Tested: sprintf (cstrName, "Power Device Tested"); break;
                case kHIDUsage_PD_AwaitingPower: sprintf (cstrName, "Power Device Awaiting Power"); break;
                case kHIDUsage_PD_CommunicationLost: sprintf (cstrName, "Power Device Communication Lost"); break;
                    
                case kHIDUsage_PD_iManufacturer: sprintf (cstrName, "Power Device Manufacturer String Index"); break;
                case kHIDUsage_PD_iProduct: sprintf (cstrName, "Power Device Product String Index"); break;
                case kHIDUsage_PD_iserialNumber: sprintf (cstrName, "Power Device Serial Number String Index"); break;
                default: sprintf (cstrName, "Power Device Usage 0x%lx", valueUsage); break;
            }
            break;
        case kHIDPage_BatterySystem:
			if (((valueUsage >= 0x0A) && (valueUsage <= 0x0F)) || ((valueUsage >= 0x1E) && (valueUsage <= 0x27)) ||
                ((valueUsage >= 0x30) && (valueUsage <= 0x3F)) || ((valueUsage >= 0x4C) && (valueUsage <= 0x5F)) ||
                ((valueUsage >= 0x6C) && (valueUsage <= 0x7F)) || ((valueUsage >= 0x90) && (valueUsage <= 0xBF)) ||
                ((valueUsage >= 0xC3) && (valueUsage <= 0xCF)) || ((valueUsage >= 0xDD) && (valueUsage <= 0xEF)) ||
                ((valueUsage >= 0xF2) && (valueUsage <= 0xFF)))
                sprintf (cstrName, "Power Device Reserved");
			else
				switch (valueUsage)
            {
                case kHIDUsage_BS_Undefined: sprintf (cstrName, "Battery System Undefined"); break;
                case kHIDUsage_BS_SMBBatteryMode: sprintf (cstrName, "SMB Mode"); break;
                case kHIDUsage_BS_SMBBatteryStatus: sprintf (cstrName, "SMB Status"); break;
                case kHIDUsage_BS_SMBAlarmWarning: sprintf (cstrName, "SMB Alarm Warning"); break;
                case kHIDUsage_BS_SMBChargerMode: sprintf (cstrName, "SMB Charger Mode"); break;
                case kHIDUsage_BS_SMBChargerStatus: sprintf (cstrName, "SMB Charger Status"); break;
                case kHIDUsage_BS_SMBChargerSpecInfo: sprintf (cstrName, "SMB Charger Extended Status"); break;
                case kHIDUsage_BS_SMBSelectorState: sprintf (cstrName, "SMB Selector State"); break;
                case kHIDUsage_BS_SMBSelectorPresets: sprintf (cstrName, "SMB Selector Presets"); break;
                case kHIDUsage_BS_SMBSelectorInfo: sprintf (cstrName, "SMB Selector Info"); break;
                case kHIDUsage_BS_OptionalMfgFunction1: sprintf (cstrName, "Battery System Optional SMB Mfg Function 1"); break;
                case kHIDUsage_BS_OptionalMfgFunction2: sprintf (cstrName, "Battery System Optional SMB Mfg Function 2"); break;
                case kHIDUsage_BS_OptionalMfgFunction3: sprintf (cstrName, "Battery System Optional SMB Mfg Function 3"); break;
                case kHIDUsage_BS_OptionalMfgFunction4: sprintf (cstrName, "Battery System Optional SMB Mfg Function 4"); break;
                case kHIDUsage_BS_OptionalMfgFunction5: sprintf (cstrName, "Battery System Optional SMB Mfg Function 5"); break;
                case kHIDUsage_BS_ConnectionToSMBus: sprintf (cstrName, "Battery System Connection To System Management Bus"); break;
                case kHIDUsage_BS_OutputConnection: sprintf (cstrName, "Battery System Output Connection Status"); break;
                case kHIDUsage_BS_ChargerConnection: sprintf (cstrName, "Battery System Charger Connection"); break;
                case kHIDUsage_BS_BatteryInsertion: sprintf (cstrName, "Battery System Battery Insertion"); break;
                case kHIDUsage_BS_Usenext: sprintf (cstrName, "Battery System Use Next"); break;
                case kHIDUsage_BS_OKToUse: sprintf (cstrName, "Battery System OK To Use"); break;
                case kHIDUsage_BS_BatterySupported: sprintf (cstrName, "Battery System Battery Supported"); break;
                case kHIDUsage_BS_SelectorRevision: sprintf (cstrName, "Battery System Selector Revision"); break;
                case kHIDUsage_BS_ChargingIndicator: sprintf (cstrName, "Battery System Charging Indicator"); break;
                case kHIDUsage_BS_ManufacturerAccess: sprintf (cstrName, "Battery System Manufacturer Access"); break;
                case kHIDUsage_BS_RemainingCapacityLimit: sprintf (cstrName, "Battery System Remaining Capacity Limit"); break;
                case kHIDUsage_BS_RemainingTimeLimit: sprintf (cstrName, "Battery System Remaining Time Limit"); break;
                case kHIDUsage_BS_AtRate: sprintf (cstrName, "Battery System At Rate..."); break;
                case kHIDUsage_BS_CapacityMode: sprintf (cstrName, "Battery System Capacity Mode"); break;
                case kHIDUsage_BS_BroadcastToCharger: sprintf (cstrName, "Battery System Broadcast To Charger"); break;
                case kHIDUsage_BS_PrimaryBattery: sprintf (cstrName, "Battery System Primary Battery"); break;
                case kHIDUsage_BS_ChargeController: sprintf (cstrName, "Battery System Charge Controller"); break;
                case kHIDUsage_BS_TerminateCharge: sprintf (cstrName, "Battery System Terminate Charge"); break;
                case kHIDUsage_BS_TerminateDischarge: sprintf (cstrName, "Battery System Terminate Discharge"); break;
                case kHIDUsage_BS_BelowRemainingCapacityLimit: sprintf (cstrName, "Battery System Below Remaining Capacity Limit"); break;
                case kHIDUsage_BS_RemainingTimeLimitExpired: sprintf (cstrName, "Battery System Remaining Time Limit Expired"); break;
                case kHIDUsage_BS_Charging: sprintf (cstrName, "Battery System Charging"); break;
                case kHIDUsage_BS_Discharging: sprintf (cstrName, "Battery System Discharging"); break;
                case kHIDUsage_BS_FullyCharged: sprintf (cstrName, "Battery System Fully Charged"); break;
                case kHIDUsage_BS_FullyDischarged: sprintf (cstrName, "Battery System Fully Discharged"); break;
                case kHIDUsage_BS_ConditioningFlag: sprintf (cstrName, "Battery System Conditioning Flag"); break;
                case kHIDUsage_BS_AtRateOK: sprintf (cstrName, "Battery System At Rate OK"); break;
                case kHIDUsage_BS_SMBErrorCode: sprintf (cstrName, "Battery System SMB Error Code"); break;
                case kHIDUsage_BS_NeedReplacement: sprintf (cstrName, "Battery System Need Replacement"); break;
                case kHIDUsage_BS_AtRateTimeToFull: sprintf (cstrName, "Battery System At Rate Time To Full"); break;
                case kHIDUsage_BS_AtRateTimeToEmpty: sprintf (cstrName, "Battery System At Rate Time To Empty"); break;
                case kHIDUsage_BS_AverageCurrent: sprintf (cstrName, "Battery System Average Current"); break;
                case kHIDUsage_BS_Maxerror: sprintf (cstrName, "Battery System Max Error"); break;
                case kHIDUsage_BS_RelativeStateOfCharge: sprintf (cstrName, "Battery System Relative State Of Charge"); break;
                case kHIDUsage_BS_AbsoluteStateOfCharge: sprintf (cstrName, "Battery System Absolute State Of Charge"); break;
                case kHIDUsage_BS_RemainingCapacity: sprintf (cstrName, "Battery System Remaining Capacity"); break;
                case kHIDUsage_BS_FullChargeCapacity: sprintf (cstrName, "Battery System Full Charge Capacity"); break;
                case kHIDUsage_BS_RunTimeToEmpty: sprintf (cstrName, "Battery System Run Time To Empty"); break;
                case kHIDUsage_BS_AverageTimeToEmpty: sprintf (cstrName, "Battery System Average Time To Empty"); break;
                case kHIDUsage_BS_AverageTimeToFull: sprintf (cstrName, "Battery System Average Time To Full"); break;
                case kHIDUsage_BS_CycleCount: sprintf (cstrName, "Battery System Cycle Count"); break;
                case kHIDUsage_BS_BattPackModelLevel: sprintf (cstrName, "Battery System Batt Pack Model Level"); break;
                case kHIDUsage_BS_InternalChargeController: sprintf (cstrName, "Battery System Internal Charge Controller"); break;
                case kHIDUsage_BS_PrimaryBatterySupport: sprintf (cstrName, "Battery System Primary Battery Support"); break;
                case kHIDUsage_BS_DesignCapacity: sprintf (cstrName, "Battery System Design Capacity"); break;
                case kHIDUsage_BS_SpecificationInfo: sprintf (cstrName, "Battery System Specification Info"); break;
                case kHIDUsage_BS_ManufacturerDate: sprintf (cstrName, "Battery System Manufacturer Date"); break;
                case kHIDUsage_BS_SerialNumber: sprintf (cstrName, "Battery System Serial Number"); break;
                case kHIDUsage_BS_iManufacturerName: sprintf (cstrName, "Battery System Manufacturer Name Index"); break;
                case kHIDUsage_BS_iDevicename: sprintf (cstrName, "Battery System Device Name Index"); break;
                case kHIDUsage_BS_iDeviceChemistry: sprintf (cstrName, "Battery System Device Chemistry Index"); break;
                case kHIDUsage_BS_ManufacturerData: sprintf (cstrName, "Battery System Manufacturer Data"); break;
                case kHIDUsage_BS_Rechargable: sprintf (cstrName, "Battery System Rechargable"); break;
                case kHIDUsage_BS_WarningCapacityLimit: sprintf (cstrName, "Battery System Warning Capacity Limit"); break;
                case kHIDUsage_BS_CapacityGranularity1: sprintf (cstrName, "Battery System Capacity Granularity 1"); break;
                case kHIDUsage_BS_CapacityGranularity2: sprintf (cstrName, "Battery System Capacity Granularity 2"); break;
                case kHIDUsage_BS_iOEMInformation: sprintf (cstrName, "Battery System OEM Information Index"); break;
                case kHIDUsage_BS_InhibitCharge: sprintf (cstrName, "Battery System Inhibit Charge"); break;
                case kHIDUsage_BS_EnablePolling: sprintf (cstrName, "Battery System Enable Polling"); break;
                case kHIDUsage_BS_ResetToZero: sprintf (cstrName, "Battery System Reset To Zero"); break;
                case kHIDUsage_BS_ACPresent: sprintf (cstrName, "Battery System AC Present"); break;
                case kHIDUsage_BS_BatteryPresent: sprintf (cstrName, "Battery System Battery Present"); break;
                case kHIDUsage_BS_PowerFail: sprintf (cstrName, "Battery System Power Fail"); break;
                case kHIDUsage_BS_AlarmInhibited: sprintf (cstrName, "Battery System Alarm Inhibited"); break;
                case kHIDUsage_BS_ThermistorUnderRange: sprintf (cstrName, "Battery System Thermistor Under Range"); break;
                case kHIDUsage_BS_ThermistorHot: sprintf (cstrName, "Battery System Thermistor Hot"); break;
                case kHIDUsage_BS_ThermistorCold: sprintf (cstrName, "Battery System Thermistor Cold"); break;
                case kHIDUsage_BS_ThermistorOverRange: sprintf (cstrName, "Battery System Thermistor Over Range"); break;
                case kHIDUsage_BS_VoltageOutOfRange: sprintf (cstrName, "Battery System Voltage Out Of Range"); break;
                case kHIDUsage_BS_CurrentOutOfRange: sprintf (cstrName, "Battery System Current Out Of Range"); break;
                case kHIDUsage_BS_CurrentNotRegulated: sprintf (cstrName, "Battery System Current Not Regulated"); break;
                case kHIDUsage_BS_VoltageNotRegulated: sprintf (cstrName, "Battery System Voltage Not Regulated"); break;
                case kHIDUsage_BS_MasterMode: sprintf (cstrName, "Battery System Master Mode"); break;
                case kHIDUsage_BS_ChargerSelectorSupport: sprintf (cstrName, "Battery System Charger Support Selector"); break;
                case kHIDUsage_BS_ChargerSpec: sprintf (cstrName, "attery System Charger Specification"); break;
                case kHIDUsage_BS_Level2: sprintf (cstrName, "Battery System Charger Level 2"); break;
                case kHIDUsage_BS_Level3: sprintf (cstrName, "Battery System Charger Level 3"); break;
                default: sprintf (cstrName, "Battery System Usage 0x%lx", valueUsage); break;
            }
            break;
        case kHIDPage_AlphanumericDisplay:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Alphanumeric Display Usage 0x%lx", valueUsage); break;
        }
            break;
		case kHIDPage_BarCodeScanner:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Bar Code Scanner Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Scale:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Scale Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_CameraControl:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Camera Control Usage 0x%lx", valueUsage); break;
        }
            break;
        case kHIDPage_Arcade:
            switch (valueUsage)
        {
            default: sprintf (cstrName, "Arcade Usage 0x%lx", valueUsage); break;
        }
            break;
        default:
			if (valueUsagePage > kHIDPage_VendorDefinedStart)
				sprintf (cstrName, "Vendor Defined Usage 0x%lx", valueUsage);
			else
				sprintf (cstrName, "Page: 0x%lx, Usage: 0x%lx", valueUsagePage, valueUsage);
            break;
    }
}

#endif
