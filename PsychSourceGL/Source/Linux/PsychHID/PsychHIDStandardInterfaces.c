/*
    PsychToolbox3/Source/Linux/PsychHID/PsychHIDStandardInterfaces.c

    PROJECTS: PsychHID only.

    PLATFORMS:  Linux.

    AUTHORS:

    mario.kleiner.de@gmail.com    mk

    HISTORY:

    27.07.2011     mk     Created.

*/

#include "PsychHIDStandardInterfaces.h"

static Display *dpy = NULL;
static Display *thread_dpy = NULL;
static int xi_opcode, event, error, major, minor;
static XIDeviceInfo *info = NULL;
static int ndevices = 0;
static int masterDevice = -1;
static XDevice* x_dev[PSYCH_HID_MAX_DEVICES];

static double* psychHIDKbQueueFirstPress[PSYCH_HID_MAX_DEVICES];
static double* psychHIDKbQueueFirstRelease[PSYCH_HID_MAX_DEVICES];
static double* psychHIDKbQueueLastPress[PSYCH_HID_MAX_DEVICES];
static double* psychHIDKbQueueLastRelease[PSYCH_HID_MAX_DEVICES];
static int*    psychHIDKbQueueScanKeys[PSYCH_HID_MAX_DEVICES];
static int     psychHIDKbQueueNumValuators[PSYCH_HID_MAX_DEVICES];
static unsigned int psychHIDKbQueueFlags[PSYCH_HID_MAX_DEVICES];
static PsychHIDEventRecord psychHIDKbQueueOldEvent[PSYCH_HID_MAX_DEVICES];
static Window  psychHIDKbQueueXWindow[PSYCH_HID_MAX_DEVICES];
static psych_bool psychHIDKbQueueActive[PSYCH_HID_MAX_DEVICES];
static psych_mutex KbQueueMutex;
static psych_condition KbQueueCondition;
static psych_bool  KbQueueThreadTerminate;
static psych_thread KbQueueThread;
static XEvent KbQueue_xevent;

static XDevice* GetXDevice(int deviceIndex)
{
    if (deviceIndex < 0 || deviceIndex >= PSYCH_HID_MAX_DEVICES) PsychErrorExitMsg(PsychError_user, "Invalid deviceIndex specified. No such device!");
    if (x_dev[deviceIndex] == NULL) x_dev[deviceIndex] = XOpenDevice(dpy, (XID) info[deviceIndex].deviceid);
    return(x_dev[deviceIndex]);
}

void PsychHIDInitializeHIDStandardInterfaces(void)
{
    int rc, i;

    // Init x_dev array:
    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) x_dev[i] = NULL;

    // Init keyboard queue arrays:
    memset(&psychHIDKbQueueFirstPress[0], 0, sizeof(psychHIDKbQueueFirstPress));
    memset(&psychHIDKbQueueFirstRelease[0], 0, sizeof(psychHIDKbQueueFirstRelease));
    memset(&psychHIDKbQueueLastPress[0], 0, sizeof(psychHIDKbQueueLastPress));
    memset(&psychHIDKbQueueLastRelease[0], 0, sizeof(psychHIDKbQueueLastRelease));
    memset(&psychHIDKbQueueActive[0], 0, sizeof(psychHIDKbQueueActive));
    memset(&psychHIDKbQueueScanKeys[0], 0, sizeof(psychHIDKbQueueScanKeys));
    memset(&psychHIDKbQueueNumValuators[0], 0, sizeof(psychHIDKbQueueNumValuators));
    memset(&psychHIDKbQueueOldEvent[0], 0, sizeof(psychHIDKbQueueOldEvent));
    memset(&psychHIDKbQueueFlags[0], 0, sizeof(psychHIDKbQueueFlags));
    memset(&psychHIDKbQueueXWindow[0], 0, sizeof(psychHIDKbQueueXWindow));

    // We must initialize XLib for multithreading-safe operations / access on first
    // call if usercode explicitely requests this via environment variable PSYCH_XINITTHREADS.
    //
    // We can only do this as opt-in, as XInitThreads() must be called as very first
    // XLib function after process startup or bad things will happen! We don't have control
    // over this wrt. Matlab or Octave (especially future Octave 3.7+ with its QT based GUI),
    // so we implemented our own locking in Screen() and don't need it in PsychHID, as PsychHID's
    // x-connection is exclusively used by PsychHID's Xinput processing thread. However, there
    // may be some cases when our own locking is insufficient, due to deficiencies in the
    // DRI2 XOrg FOSS Mesa graphics driver stack, so some users may want to opt-into use
    // XLib's threading protection as a work-around if they can guarantee Octave or Matlab
    // hasn't called any XLib calls already during its running session:
    if (getenv("PSYCH_XINITTHREADS")) XInitThreads();

    // Open our own private X-Display connection for HID handling:
    dpy = XOpenDisplay(NULL);
    if (!dpy) {
        PsychErrorExitMsg(PsychError_system, "PsychHID: FATAL ERROR: Couldn't open default X11 display connection to X-Server! Game over!");
    }

    // XInputExtension supported? If so do basic init:
    if (!XQueryExtension(dpy, "XInputExtension", &xi_opcode, &event, &error)) {
        printf("PsychHID: ERROR: XINPUT extension unsupported!\n");
        goto out;
    }

    // XInput V2.2 supported?
    major = 2;
    minor = 2;
    rc = XIQueryVersion(dpy, &major, &minor);
    if (rc == BadRequest) {
        printf("PsychHID: ERROR: No XI2 support. Server supports version %d.%d only.\n", major, minor);
        goto out;
    } else if (rc != Success) {
        printf("PsychHID: ERROR: Internal Error! This is a bug in Xlib.\n");
        goto out;
    }

    // Open our own 2nd private X-Display connection for HID handling. This one is for the
    // async KbQueue processing thread:
    thread_dpy = XOpenDisplay(NULL);
    if (!thread_dpy) {
        PsychErrorExitMsg(PsychError_system, "PsychHID: FATAL ERROR: Couldn't open default X11 display connection to X-Server! Game over!");
    }

    // XInputExtension supported? If so do basic init:
    if (!XQueryExtension(thread_dpy, "XInputExtension", &xi_opcode, &event, &error)) {
        printf("PsychHID: ERROR: XINPUT extension unsupported!\n");
        goto out;
    }

    // XInput V2.2 supported?
    major = 2;
    minor = 2;
    rc = XIQueryVersion(thread_dpy, &major, &minor);
    if (rc == BadRequest) {
        printf("PsychHID: ERROR: No XI2 support. Server supports version %d.%d only.\n", major, minor);
        goto out;
    } else if (rc != Success) {
        printf("PsychHID: ERROR: Internal Error! This is a bug in Xlib.\n");
        goto out;
    }

    // Enumerate all XI2 input devices:
    info = XIQueryDevice(dpy, XIAllDevices, &ndevices);

    masterDevice = -1;
    for(i = 0; i < ndevices; i++) {
        if (info[i].use == XIMasterPointer) {
            masterDevice = i;
            break;
        }
    }
    if (masterDevice == -1) printf("PsychHID: WARNING! No master keyboard/pointer found! This will end badly...\n");

    // Switch X-Server connection into synchronous mode: We need this to get
    // a higher timing precision.
    XSynchronize(dpy, TRUE);

    // Create keyboard queue mutex for later use:
    KbQueueThreadTerminate = FALSE;
    PsychInitMutex(&KbQueueMutex);
    PsychInitCondition(&KbQueueCondition, NULL);

    return;

out:
    if (dpy) XCloseDisplay(dpy);
    dpy = NULL;

    if (thread_dpy) XCloseDisplay(thread_dpy);
    thread_dpy = NULL;

    PsychErrorExitMsg(PsychError_system, "PsychHID: FATAL ERROR: X Input extension version 2.0 or later not available! Game over!");
}

void PsychHIDShutdownHIDStandardInterfaces(void)
{
    int i;

    // Close all devices registered in x_dev array:
    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
        if (x_dev[i]) XCloseDevice(dpy, x_dev[i]);
        x_dev[i] = NULL;
    }

    // Release all keyboard queues:
    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
        if (psychHIDKbQueueFirstPress[i]) {
            PsychHIDOSKbQueueRelease(i);
        }
    }

    // Release keyboard queue mutex:
    PsychDestroyMutex(&KbQueueMutex);
    PsychDestroyCondition(&KbQueueCondition);
    KbQueueThreadTerminate = FALSE;

    // Release list of enumerated input devices:
    XIFreeDeviceInfo(info);
    info = NULL;

    // Close our dedicated x-display connection and we are done:
    if (dpy) XCloseDisplay(dpy);
    dpy = NULL;

    if (thread_dpy) XCloseDisplay(thread_dpy);
    thread_dpy = NULL;

    return;
}

// Return -1 if not a touch device. Return >= 0 for touch device.
// Returns 0 for "unknown number of touch points - device didn't tell us".
// Otherwise returns > 0 for number of maximally supported touch points.
// *type is 0 for dependent touch devices like trackpads, 1 for real touch screens.
int PsychHIDIsTouchDevice(int deviceIndex, int* type)
{
    int j, count = 0, num_touches = 0;
    XIDeviceInfo *dev = &info[deviceIndex];

    if (type)
        *type = -1;

    // XInput 2.2+ supported? Otherwise no touch device support.
    if (minor >= 2) {
        for (j = 0; j < dev->num_classes; j++) {
            XIAnyClassInfo *class = dev->classes[j];
            XITouchClassInfo *t = (XITouchClassInfo*) class;

            if (class->type != XITouchClass)
                continue;

            count++;
            num_touches += t->num_touches;

            if (type)
                *type = (t->mode == XIDirectTouch) ? 1 : 0;
        }
    }

    return((count > 0) ? num_touches : -1);
}

PsychError PsychHIDEnumerateHIDInputDevices(int deviceClass)
{
    const char *deviceFieldNames[]={"usagePageValue", "usageValue", "usageName", "index", "transport", "vendorID", "productID", "version",
                                    "manufacturer", "product", "serialNumber", "locationID", "interfaceID", "totalElements", "features", "inputs",
                                    "outputs", "collections", "axes", "buttons", "hats", "sliders", "dials", "wheels", "touchDeviceType", "maxTouchpoints"};
    int numDeviceStructElements, numDeviceStructFieldNames=26, deviceIndex;
    PsychGenericScriptType *deviceStruct;
    XIDeviceInfo *dev;
    int i, j;
    int numKeys, numAxis, touchType;
    char *type = "";

    // Preparse: Count matching devices for deviceClass
    numDeviceStructElements = 0;
    for(i = 0; i < ndevices; i++) {
        dev = &info[i];
        if ((int) (dev->use) == deviceClass) numDeviceStructElements++;
    }

    // Alloc struct array of sufficient size:
    PsychAllocOutStructArray(1, FALSE, numDeviceStructElements, numDeviceStructFieldNames, deviceFieldNames, &deviceStruct);
    deviceIndex = 0;

    // Return info:
    for(i = 0; i < ndevices; i++) {
        // Check i'th device:
        dev = &info[i];

        // Skip if non matching class:
        if ((int) (dev->use) != deviceClass) continue;

        switch(dev->use) {
            case XIMasterPointer: type = "master pointer"; break;
            case XIMasterKeyboard: type = "master keyboard"; break;
            case XISlavePointer: type = "slave pointer"; break;
            case XISlaveKeyboard: type = "slave keyboard"; break;
            case XIFloatingSlave: type = "floating slave"; break;
        }

        // Usagepage is 1 for "Desktop usage page":
        PsychSetStructArrayDoubleElement("usagePageValue", deviceIndex,  (double) 1, deviceStruct);

        if (dev->use == XIMasterKeyboard || dev->use == XISlaveKeyboard) {
            // Usage 6 is for keyboard:
            PsychSetStructArrayDoubleElement("usageValue", deviceIndex, (double) 6, deviceStruct);
        }
        else if (dev->use != XIFloatingSlave) {
            // Usage is 2 for mouse:
            PsychSetStructArrayDoubleElement("usageValue", deviceIndex, (double) 2, deviceStruct);
        } else {
            // Assign meaningless Usage of zero for floating slave devices:
            PsychSetStructArrayDoubleElement("usageValue", deviceIndex, (double) 0, deviceStruct);
        }

        PsychSetStructArrayStringElement("usageName",  deviceIndex,  type, deviceStruct);
        PsychSetStructArrayDoubleElement("index",  deviceIndex,  (double) i, deviceStruct);
        PsychSetStructArrayStringElement("transport",  deviceIndex,  (dev->enabled) ? "enabled" : "disabled", deviceStruct);
        PsychSetStructArrayStringElement("product",  deviceIndex,  dev->name, deviceStruct);
        PsychSetStructArrayDoubleElement("locationID",  deviceIndex,  (double) dev->attachment, deviceStruct);
        PsychSetStructArrayDoubleElement("interfaceID",  deviceIndex,  (double) dev->deviceid, deviceStruct);

        //PsychSetStructArrayDoubleElement("vendorID",  deviceIndex,  (double)currentDevice->vendorID,  deviceStruct);
        //PsychSetStructArrayDoubleElement("productID",  deviceIndex,  (double)currentDevice->productID,  deviceStruct);
        //PsychSetStructArrayDoubleElement("version",  deviceIndex,  (double)currentDevice->version,  deviceStruct);
        //PsychSetStructArrayStringElement("manufacturer", deviceIndex,  currentDevice->manufacturer,   deviceStruct);
        //PsychSetStructArrayStringElement("serialNumber", deviceIndex,  currentDevice->serial,    deviceStruct);

        numKeys = numAxis = 0;
        for (j = 0; j < dev->num_classes; j++) {
            if (dev->classes[j]->type == XIKeyClass) numKeys += (int) (((XIKeyClassInfo*) dev->classes[j])->num_keycodes);
            if (dev->classes[j]->type == XIButtonClass) numKeys += (int) (((XIButtonClassInfo*) dev->classes[j])->num_buttons);
            if (dev->classes[j]->type == XIValuatorClass) numAxis++;
        }

        PsychSetStructArrayDoubleElement("totalElements", deviceIndex,  (double) numKeys + numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("features",  deviceIndex,  (double) dev->num_classes, deviceStruct);
        PsychSetStructArrayDoubleElement("inputs",  deviceIndex,  (double) numKeys + numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("outputs",  deviceIndex,  (double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("collections",         deviceIndex,  (double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("axes",  deviceIndex,  (double) numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("buttons",  deviceIndex,  (double) numKeys, deviceStruct);
        PsychSetStructArrayDoubleElement("hats",  deviceIndex,  (double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("sliders",  deviceIndex,  (double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("dials",  deviceIndex,  (double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("wheels",  deviceIndex,  (double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("maxTouchpoints",  deviceIndex, (double) PsychHIDIsTouchDevice(i, NULL), deviceStruct);
        PsychHIDIsTouchDevice(i, &touchType);
        PsychSetStructArrayDoubleElement("touchDeviceType",  deviceIndex, (double) touchType, deviceStruct);

        deviceIndex++;
    }

    return(PsychError_none);
}

PsychError PsychHIDOSKbCheck(int deviceIndex, double* scanList)
{
    double* buttonStates;
    unsigned char keys_return[32];
    int keysdown;
    double timestamp;
    int i, j;
    psych_bool isButtons = FALSE;

    memset(keys_return, 0, sizeof(keys_return));

    // Map "default" deviceIndex to legacy "Core protocol" method of querying keyboard
    // state. This will give us whatever X has setup as default keyboard:
    if (deviceIndex == INT_MAX) {
        // Request current keyboard state of default keyboard from X-Server:
        XQueryKeymap(dpy, (char*) keys_return);
    } else if (deviceIndex < 0 || deviceIndex >= ndevices) {
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard deviceIndex specified. No such device!");
    } else if (info[deviceIndex].use == XIMasterKeyboard) {
        // Master keyboard:

        // Query current client pointer assignment, then switch it to
        // associated master pointer for the master keyboard we want
        // to query. This way, all future queries will query our requested
        // master keyboard:
        j = -1;
        if (!XIGetClientPointer(dpy, None, &j) || (j != info[deviceIndex].attachment)) XISetClientPointer(dpy, None, info[deviceIndex].attachment);

        // Request current keyboard state from X-Server:
        XQueryKeymap(dpy, (char*) keys_return);

        // Reset master pointer/keyboard assignment to pre-query state:
        if ((j > 0) && (j != info[deviceIndex].attachment)) XISetClientPointer(dpy, None, j);
    } else {
        // Non-Default deviceIndex: Want to query specific slave keyboard.
        if (info[deviceIndex].use == XIMasterPointer) PsychErrorExitMsg(PsychError_user, "Invalid deviceIndex specified! Cannot query master mouse pointers as keyboards.");

        // Open connection to non-master-keyboard device:
        XDevice* mydev = GetXDevice(deviceIndex);

        // Query its current state:
        XDeviceState* state = XQueryDeviceState(dpy, mydev);
        XInputClass* data = state->data;

        // printf("Dummy = %i , NClasses = %i\n", dummy1, state->num_classes);

        // Find state structure with key status info:
        for (i = 0; i < state->num_classes; i++) {
            // printf("Class %i: Type %i - %i\n", i, (int) data->class, (int) data->length);
            if (data->class == KeyClass) {
                // printf("NumKeys %i\n", ((XKeyState*) data)->num_keys);

                // Copy 32 Byte keystate vector into key_return. Each bit encodes for one key:
                memcpy(&keys_return[0], &(((XKeyState*) data)->keys[0]), sizeof(keys_return));
                isButtons = FALSE;
            }

            // Also handle devices with buttons as if they are keyboards, e.g., mouse, joystick...
            if (data->class == ButtonClass) {
                // printf("NumButtons %i\n", ((XButtonState*) data)->num_buttons);

                // Copy 32 Byte buttonstate vector into key_return. Each bit encodes for one button:
                memcpy(&keys_return[0], &(((XButtonState*) data)->buttons[0]), sizeof(keys_return));
                isButtons = TRUE;
            }

            // Advance to next entry:
            data = (XInputClass*) (((void*) data) + ((size_t) data->length));
        }

        XFreeDeviceState(state);
    }

    // Done with query. We have keyboard state in keys_return[] now.

    // Request current time of query:
    PsychGetAdjustedPrecisionTimerSeconds(&timestamp);

    // Reset overall key state to "none pressed":
    keysdown = 0;

    // Any key down?
    for (i = 0; i < 32; i++) keysdown+=(unsigned int) keys_return[i];

    // Copy out overall keystate:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (keysdown > 0) ? 1 : 0);

    // Copy out timestamp:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, timestamp);

    // Copy keyboard state:
    PsychAllocOutDoubleMatArg(3, kPsychArgOptional, 1, 256, 1, &buttonStates);

    // Map 32 times 8 bitvector to 256 element return vector:
    for(i = 0; i < 32; i++) {
        for(j = 0; j < 8; j++) {
            // This button or key down?
            buttonStates[i*8 + j] = (keys_return[i] & (1<<j)) ? 1 : 0;

            // Is this a button, e.g., on mouse/joystick, instead of a key on keyboard?
            if (isButtons) {
                // All buttons are shifted to + 1 index, so we need to shift back by
                // 1 position. Do this in a pathetic way to save brain-cycles:
                if (i*8 + j > 0) buttonStates[i*8 + j - 1] = buttonStates[i*8 + j];
                if (i*8 + j >= 255) buttonStates[i*8 + j] = 0;
            }

            // Apply scanList mask, if any provided:
            if (scanList && (scanList[i*8 + j] <= 0)) buttonStates[i*8 + j] = 0;
        }
    }

    return(PsychError_none);
}

PsychError PsychHIDOSGamePadAxisQuery(int deviceIndex, int axisId, double* min, double* max, double* val, char* axisLabel)
{
    XIDeviceInfo *dev = NULL;
    int i, j, dummy1;

    dev = XIQueryDevice(dpy, info[deviceIndex].deviceid, &dummy1);

    printf("Dummy = %i , NClasses = %i\n", dummy1, dev->num_classes);
    for (i = 0; i < dev->num_classes; i++) {
        printf("Class %i: Type %i\n", i, (int) dev->classes[i]->type);
        if (dev->classes[i]->type == XIKeyClass) printf("NumKeys %i\n", ((XIKeyClassInfo*) dev->classes[i])->num_keycodes);
        if (dev->classes[i]->type == XIButtonClass) printf("NumButtons %i\n", ((XIButtonClassInfo*) dev->classes[i])->num_buttons);
        if (dev->classes[i]->type == XIValuatorClass) {
            printf("Value %f\n", (float) ((XIValuatorClassInfo*) dev->classes[i])->value);
            if (val) *val = (double) ((XIValuatorClassInfo*) dev->classes[i])->value;
        }
    }

    XIFreeDeviceInfo(dev);

    // Open connection to slave keyboard device:
    XDevice* mydev = GetXDevice(deviceIndex);

    // Query its current state:
    XDeviceState* state = XQueryDeviceState(dpy, mydev);

    printf("NClasses = %i\n", state->num_classes);

    // Find state structure with key status info:
    XInputClass* data = state->data;
    for (i = 0; i < state->num_classes; i++) {
        printf("Class = %i\n", (int) data->class);
        if (data->class == ValuatorClass) {
            XValuatorState* valuator = (XValuatorState*) data;
            printf("NumAxis %i [Mode=%s]\n", valuator->num_valuators, (valuator->mode == Absolute) ? "Absolute" : "Relative");
            for (j = 0 ; j < valuator->num_valuators; j++) printf("Axis %i: %i\n", j, valuator->valuators[j]);
        }

        data = (XInputClass*) (((void*) data) + ((size_t) data->length));
    }

    XFreeDeviceState(state);

    int nSamples, mode_return, axis_count_return;
    double* outSamples;
    double tSample;

    // Query motion history of this device:
    XDeviceTimeCoord* samples = XGetDeviceMotionEvents(dpy, mydev, 0, CurrentTime, &nSamples, &mode_return, &axis_count_return);
    printf("[%p] n = %i , mr = %s , ac = %i\n", samples, nSamples, (mode_return == Absolute) ? "Absolute" : "Relative" , axis_count_return);

    // Return it as 1st argument to userspace if requested:
    if (PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 2 + axis_count_return, nSamples, 1, &outSamples)) {
        for (i = 0; i < nSamples; i++) {
            // Sampleindex in row 1:
            *(outSamples++) = (double) i;

            // Sampletime in row 2: Need to convert msecs to sec and
            // then map from CLOCK_MONOTONIC to our CLOCK_REALTIME
            // GetSecs() timebase:
            tSample = (double) samples[i].time / 1000.0;
            *(outSamples++) = PsychOSMonotonicToRefTime(tSample);

            // Axis motion samples in successive rows:
            for (j = 0; j < axis_count_return; j++) {
                *(outSamples++) = (double) samples[i].data[j];
            }
            // Next sample in next column:
        }
    }

    // Release returned sample motion buffer:
    XFreeDeviceMotionEvents(samples);

    return(PsychError_none);
}

// This is the event dequeue & process function which updates
// Keyboard queue state. It can be called with 'blockingSinglepass'
// set to TRUE to process exactly one event, if called from the
// background keyboard queue processing thread. Alternatively it
// can be called synchronously from KbQueueCheck with a setting of FALSE
// to iterate over all available events and process them instantaneously:
void KbQueueProcessEvents(psych_bool blockingSinglepass)
{
    PsychHIDEventRecord evt;
    XKeyPressedEvent key;
    XIDeviceEvent* event;
    XIRawEvent* rawevent;
    psych_bool valid;
    double tnow;
    int i, j, index, deviceid, numValuators;
    unsigned int screen_width, screen_height;
    char asciiChar;

    while (1) {
        XGenericEventCookie *cookie = &KbQueue_xevent.xcookie;

        // Single pass or multi-pass?
        if (blockingSinglepass) {
            // Wait until at least one event available and dequeue it:
            XNextEvent(thread_dpy, &KbQueue_xevent);
        } else {
            // Check if event available, dequeue it, if so. Abort
            // processing if no new event available, aka queue empty:
            if (!XCheckTypedEvent(thread_dpy, GenericEvent, &KbQueue_xevent)) break;
        }

        // Take timestamp:
        PsychGetAdjustedPrecisionTimerSeconds(&tnow);

        // Clear ringbuffer event:
        memset(&evt, 0 , sizeof(evt));

        // Set cooked character to "undefined" as a starter: It will only get
        // assigned something else if keypress/release events from real keyboards
        // or keypads are processed:
        evt.cookedEventCode = -1;

        // Is this an event we're interested in?
        if ((cookie->type == GenericEvent) && (cookie->extension == xi_opcode)) {
            // Yes. Process it:
            if (XGetEventData(thread_dpy, cookie)) {
                // Process it:

                // printf("Event type %d received\n", cookie->evtype);
                if (cookie->evtype == XI_RawButtonPress || cookie->evtype == XI_RawButtonRelease || cookie->evtype == XI_RawMotion ||
                    cookie->evtype == XI_RawTouchBegin || cookie->evtype == XI_RawTouchEnd || cookie->evtype == XI_RawTouchUpdate) {
                    // Raw device event for mice and similar devices:
                    rawevent = (XIRawEvent*) cookie->data;
                    event = NULL;
                    valid = TRUE; // Always true for raw devices like mice etc., unless queue flag 1 is set, see below.
                    index = rawevent->detail;
                    deviceid = rawevent->deviceid;
                }
                else {
                    // Regular device event:
                    Window rootRet;
                    unsigned int depth_return, border_width_return;
                    int x, y;

                    event = (XIDeviceEvent*) cookie->data;
                    rawevent = NULL;
                    valid = !(event->flags & XIKeyRepeat);
                    index = event->detail;
                    deviceid = event->deviceid;

                    // Get width x height of associated root window, aka screen size for touch coordinate remapping:
                    XGetGeometry(thread_dpy, event->root, &rootRet, &x, &y, &screen_width, &screen_height, &border_width_return, &depth_return);
                }

                if (event && (cookie->evtype != XI_TouchOwnership)) {
                    // Process device button state:
                    for (j = 0; (j < event->buttons.mask_len * 8) && (j < 32); j++) {
                        // j'th device button pressed?
                        if (XIMaskIsSet(event->buttons.mask, j)) {
                            // Yes: Set status to "pressed" as at least one device button is pressed:
                            // Note this gets overriden for some event types below, like touch and
                            // keyboard. It only really affects mouse/gamepad/joystick movements.
                            evt.status |= (1 << 0);

                            // Mark button pressed:
                            evt.buttonStates |= (1 << j);
                        }
                    }
                }

                // Map Xinput device id to PTB 'deviceIndex' aka the proper keyboard queue:
                for (i = 0; i < ndevices; i++) if (deviceid == info[i].deviceid) break;

                // Special handling for synthetic key repeat flags required?
                if ((i < ndevices) && (psychHIDKbQueueFlags[i] & 0x3)) {
                    // Yes: Filter out key repeat on raw events (from pointer devices) as well?
                    // This will suppress scroll events (mouse wheel, PowerMate knob etc.):
                    if (rawevent && (psychHIDKbQueueFlags[i] & 0x1))
                        valid = !(rawevent->flags & XIKeyRepeat);

                    // Always accept key repeats:
                    if (psychHIDKbQueueFlags[i] & 0x2)
                        valid = TRUE;
                }

                // We're only interested in key press and release events, and only in
                // real ones, not XServer generated synthetic key auto-repeat events.
                // Also only for a device that we've registered:
                if ((i < ndevices) && valid &&
                    ((cookie->evtype == XI_KeyPress) || (cookie->evtype == XI_KeyRelease) ||
                     (cookie->evtype == XI_ButtonPress) || (cookie->evtype == XI_ButtonRelease) ||
                     (cookie->evtype == XI_RawButtonPress) || (cookie->evtype == XI_RawButtonRelease))) {

                    // If this is a button event from a mouse/joystick etc. instead of a key event from a keyboard/keypad, then
                    // all button indices are shifted by an offset of 1 for some weird reason. Decrement index by one to compensate
                    // for this. [Tested on Ubuntu 10.10 and 11.10 with two mice and 1 joystick]
                    if ((index > 0) &&
                        ((cookie->evtype == XI_ButtonPress) || (cookie->evtype == XI_ButtonRelease) ||
                         (cookie->evtype == XI_RawButtonPress) || (cookie->evtype == XI_RawButtonRelease))) {
                        index--;
                    }

                    // Key release on keyboard maps to character code 0.
                    if (cookie->evtype == XI_KeyRelease) evt.cookedEventCode = 0;

                    // Key press on a real keyboard needs to be mapped to character ascii code, if possible:
                    if (cookie->evtype == XI_KeyPress) {
                        // Assign info from our XIDeviceEvent to a standard XKeyPressedEvent which
                        // XLookupString() can actually understand:
                        key.type         = KeyPress;
                        key.root         = event->root;
                        key.window       = event->event;
                        key.subwindow    = event->child;
                        key.time         = event->time;
                        key.x            = event->event_x;
                        key.y            = event->event_y;
                        key.x_root       = event->root_x;
                        key.y_root       = event->root_y;
                        key.same_screen  = True;
                        key.send_event   = False;
                        key.serial       = event->serial;
                        key.display      = thread_dpy;

                        key.keycode      = index;
                        key.state        = event->mods.effective;

                        if (1 == XLookupString(((XKeyEvent*) (&key)), &asciiChar, 1, NULL, NULL)) {
                            // Mapped: Assign it.
                            evt.cookedEventCode = (int) asciiChar;
                        }
                        else {
                            // Not mappable:
                            evt.cookedEventCode = 0;
                        }

                        // Was this a CTRL + C interrupt request? This is mapped to ASCII control character 3 "ETX".
                        if (evt.cookedEventCode == 3) {
                            // Yes: Tell ConsoleInputHelper() about it, to reenable keystroke
                            // character dispatch in the terminal. This will undo a potential ListenChar(2)
                            // op if GetChar() et al. are used with the help of this keyboard queue from
                            // within Octave or Matlab in -nojvm mode.
                            ConsoleInputHelper(-1);
                        }

                        if (evt.cookedEventCode >= 0) {
                            // Tell ConsoleInputHelper() about the character:
                            ConsoleInputHelper(evt.cookedEventCode);
                        }
                    }

                    // Need the lock from here on:
                    PsychLockMutex(&KbQueueMutex);

                    // This keyboard queue created and started? Interested in this
                    // keycode?
                    if (psychHIDKbQueueActive[i] && (psychHIDKbQueueScanKeys[i][index] != 0)) {
                        // Yes: The queue wants to receive info about this key event.

                        // Press or release?
                        if ((cookie->evtype == XI_KeyPress) || (cookie->evtype == XI_ButtonPress) || (cookie->evtype == XI_RawButtonPress)) {
                            // Enqueue key press. Always in the "last press" array, because any
                            // press at this time is the best candidate for the last press.
                            // Only enqeue in "first press" if there wasn't any registered before,
                            // ie., the slot is so far empty:
                            if (psychHIDKbQueueFirstPress[i][index] == 0) psychHIDKbQueueFirstPress[i][index] = tnow;
                            psychHIDKbQueueLastPress[i][index] = tnow;
                            evt.status |= (1 << 0);
                        } else {
                            // Enqueue key release. See logic above:
                            if (psychHIDKbQueueFirstRelease[i][index] == 0) psychHIDKbQueueFirstRelease[i][index] = tnow;
                            psychHIDKbQueueLastRelease[i][index] = tnow;
                            evt.status &= ~(1 << 0);
                        }

                        // Update event buffer:
                        evt.timestamp = tnow;
                        evt.rawEventCode = index + 1;

                        // Also provide (x,y) position of the pointing device at time of event:
                        // Note: For raw events no such absolute position info is available.
                        if (event) {
                            evt.X = event->event_x;
                            evt.Y = event->event_y;
                            evt.normX = evt.X / screen_width;
                            evt.normY = evt.Y / screen_height;
                        }
                        else {
                            evt.X = evt.Y = evt.normX = evt.normY = 0;
                        }

                        PsychHIDAddEventToEventBuffer(i, &evt);

                        // Tell waiting userspace (under KbQueueMutex protection for better scheduling) something interesting has changed:
                        PsychSignalCondition(&KbQueueCondition);
                    }

                    // Done with shared data access:
                    PsychUnlockMutex(&KbQueueMutex);
                }
                else {
                    PsychLockMutex(&KbQueueMutex);

                    if ((i < ndevices) && valid && psychHIDKbQueueActive[i] && (psychHIDKbQueueNumValuators[i] > 0)) {
                        // Handling of motion events / valuator change events / touch events. These only go into the event buffer,
                        // not legacy KbQueue arrays, as those don't make sense here:
                        numValuators = psychHIDKbQueueNumValuators[i];

                        evt.cookedEventCode = -1;
                        evt.timestamp = tnow;
                        evt.rawEventCode = 0;

                        // Mouse/Touchpad/Trackpad/Trackpoint/Trackball/Joystick/Pointing device
                        // change of position or axis state:
                        if ((cookie->evtype == XI_Motion || cookie->evtype == XI_RawMotion) && (numValuators >= 2)) {
                            // Type is a pointer or axis motion event:
                            evt.type = 1;

                            // Motion event:
                            evt.status |= (1 << 1);

                            // Store absolute position:
                            // Note: For raw events no such absolute position info is available.
                            if (event) {
                                evt.X = event->event_x;
                                evt.Y = event->event_y;
                                evt.normX = evt.X / screen_width;
                                evt.normY = evt.Y / screen_height;
                            }
                            else {
                                evt.X = evt.Y = evt.normX = evt.normY = 0;
                            }

                            // More than 2 standard valuators for (x,y) position or deflection?
                            if (numValuators >= 2) {
                                if (event) {
                                    // Store up to numValuators valuator values:
                                    double *valuator = event->valuators.values;
                                    for (j = 0; (j < event->valuators.mask_len * 8) && (j < numValuators); j++) {
                                        evt.numValuators++;
                                        // Updated valuator value?
                                        if (XIMaskIsSet(event->valuators.mask, j)) {
                                            // Yes: Assign.
                                            evt.valuators[j] = (float) *valuator;
                                            valuator++;
                                        }
                                        else {
                                            // No: Assign old value from last pass:
                                            evt.valuators[j] = (float) psychHIDKbQueueOldEvent[i].valuators[j];
                                        }
                                    }
                                }
                                else {
                                    // Raw motion event:
                                    // Store up to numValuators raw values:
                                    double *raw_values = rawevent->raw_values;
                                    for (j = 0; (j < rawevent->valuators.mask_len * 8) && (j < numValuators); j++) {
                                        evt.numValuators++;
                                        // Updated valuator value?
                                        if (XIMaskIsSet(rawevent->valuators.mask, j)) {
                                            // Yes: Assign.
                                            evt.valuators[j] = (float) *raw_values;
                                            raw_values++;
                                        }
                                        else {
                                            // No: Assign old value from last pass:
                                            evt.valuators[j] = (float) psychHIDKbQueueOldEvent[i].valuators[j];
                                        }
                                    }

                                    // Derive somewhat X, Y values from first two valuators:
                                    evt.X = evt.valuators[0];
                                    evt.Y = evt.valuators[1];
                                }

                                // Keep track of last event for above valuator updating:
                                memcpy(&psychHIDKbQueueOldEvent[i], &evt, sizeof(evt));
                            }
                        }

                        if ((cookie->evtype == XI_TouchBegin || cookie->evtype == XI_TouchEnd || cookie->evtype == XI_TouchUpdate || cookie->evtype == XI_TouchOwnership ||
                             cookie->evtype == XI_RawTouchBegin || cookie->evtype == XI_RawTouchEnd || cookie->evtype == XI_RawTouchUpdate) && (numValuators >= 4)) {
                            XIDeviceInfo *dev = &info[i];

                            // Touch point id, 32-Bit integer - unique until 32 bit wraparound ;) :
                            evt.rawEventCode = (event) ? event->detail : rawevent->detail;

                            // Fetch most recent touch record in the series for this touch point:
                            // oldevt == NULL if none yet exists, or none exists anymore due to some buffer wraparound:
                            PsychHIDEventRecord *oldevt = PsychHIDLastTouchEventFromEventBuffer(i, evt.rawEventCode);

                            // Everything of interest is in the valuators:
                            if (cookie->evtype != XI_TouchOwnership) {
                                // Store up to numValuators valuator values:
                                if (event) {
                                    double *valuator = event->valuators.values;
                                    for (j = 0; (j < event->valuators.mask_len * 8) && (j < numValuators); j++) {
                                        evt.numValuators++;

                                        // Updated valuator value?
                                        if (XIMaskIsSet(event->valuators.mask, j)) {
                                            // Yes: Assign.
                                            evt.valuators[j] = (float) *valuator;
                                            valuator++;
                                        }
                                        else {
                                            // No: Assign old value from last pass:
                                            evt.valuators[j] = (float) ((oldevt) ? oldevt->valuators[j] : 0.0);
                                        }
                                    }
                                }
                                else {
                                    double *valuator = rawevent->raw_values;
                                    for (j = 0; (j < rawevent->valuators.mask_len * 8) && (j < numValuators); j++) {
                                        evt.numValuators++;

                                        // Updated valuator value?
                                        // TODO: Workaround for what might be a X-Server bug as of 1.19.6.
                                        // XIMaskIsSet() reports true also for axis/valuators that don't have
                                        // actual updated values! The corresponding *valuator value is always
                                        // zero if the axis didn't update in this event!
                                        // Work around this by only updating for non-zero *valuator values.
                                        // Use cached value from previous touch point update otherwise.
                                        // Downside is obviously that we can't report any zero values for a valuator.
                                        //
                                        // This is fine for (x,y) axis, as a user hitting exactly x==0 or y==0
                                        // ie., the left/upper touchscreen edges, is unlikely. We might get away with it
                                        // as well for pressure or area valuators, as they only become zero if
                                        // a touch ends. Valuators for touch orientation or such, where 0 degrees is
                                        // a valid reported angle, will not work properly though. This needs some fix
                                        // in the upstream X-Server, but for the moment this is the best we can do:
                                        // Disabled: if (XIMaskIsSet(rawevent->valuators.mask, j)) {
                                        if (XIMaskIsSet(rawevent->valuators.mask, j) && *valuator) {
                                            // Yes: Assign.
                                            // printf("%i: Valuator[%i]: %f\n", evt.rawEventCode, j, *valuator);
                                            evt.valuators[j] = (float) *valuator;
                                            valuator++;
                                        }
                                        else {
                                            // No: Assign old value from last pass:
                                            // printf("%i: Valuator[%i]: Gap\n", evt.rawEventCode, j);
                                            evt.valuators[j] = (float) ((oldevt) ? oldevt->valuators[j] : 0.0);
                                        }
                                    }
                                }
                            }

                            // Normalize x and y positions to 0.0 - 1.0 range (at least for absolute touch devices).
                            // Also remap to X-Screen space:
                            for (j = 0; j < dev->num_classes; j++) {
                                XIAnyClassInfo *class = dev->classes[j];
                                XIValuatorClassInfo *v = (XIValuatorClassInfo*) class;
                                double r;

                                if (class->type != XIValuatorClass)
                                    continue;

                                r = v->max - v->min;
                                if (r > 0) {
                                    if (v->number == 0)
                                        evt.normX = (float) ((evt.valuators[0] - v->min) / r);

                                    if (v->number == 1)
                                        evt.normY = (float) ((evt.valuators[1] - v->min) / r);
                                }

                                // Remap to X-Screen / RootWindow width x height:
                                evt.X = evt.normX * screen_width;
                                evt.Y = evt.normY * screen_height;
                            }

                            switch (cookie->evtype) {
                                case XI_TouchBegin:
                                case XI_RawTouchBegin:
                                    // Touch begins - Finger pressed, integrity bit cleared:
                                    evt.status = 1;
                                    evt.type = 2;
                                    break;

                                case XI_TouchUpdate:
                                case XI_RawTouchUpdate:
                                    // Finger motion event, with finger pressed onto the touch surface:
                                    evt.status |= 1 + 2;
                                    evt.type = 3;

                                    // If we get signalled that this is the end of the physical touch sequence,
                                    // but somebody else has a grab on it and hasn't decided to accept/reject
                                    // yet, then we know we've seen the full touch sequence. We may or may not
                                    // get a XI_TouchOwnership event before the XI_TouchEnd, but we know we are
                                    // good, so simply set the integrity bit:
                                    //
                                    // See http://who-t.blogspot.de/2012/01/multitouch-in-x-touch-grab-handling.html
                                    if ((event && event->flags & XITouchPendingEnd) || (rawevent && rawevent->flags & XITouchPendingEnd)) {
                                        // Set integrity bit for this touch point:
                                        evt.status |= (1 << 31);
                                        // printf("TouchPendingEnd!!\n");
                                    }
                                    break;

                                case XI_TouchEnd:
                                case XI_RawTouchEnd:
                                    // Touch ends - Finger released and no more motion:
                                    evt.status &= ~(1 + 2);
                                    evt.type = 4;

                                    // By this time we must have the integrity bit set, or we lost data.
                                    // TODO: No actually! This produces lots of false positives for multi-touch sequences,
                                    // because usually only the 1st (primary) touch in a sequence will receive the XI_TouchOwnership event
                                    // and so gets the integrity bit set on its touch id. The same is true for XITouchPendingEnd marking
                                    // on only one XI_TouchUpdate event for one touchpoint of a sequence. If multi-touches other than those
                                    // marked ones are released and ended, then this will signal a false positive touch sequence data loss.
                                    // Additionally, if raw touch events are used, we don't seem to receive XI_TouchOwnership events at all.
                                    //
                                    // We either need to do without sequence loss detection or completely rethink the approach.
                                    // Leave it here for now, but disable loss reporting:
                                    //if (!oldevt || !(oldevt->status & (1 << 31))) {
                                    if (FALSE) {
                                        // Nope, we lost touch data. Did we already send a fail event for this touch queue?
                                        // If not, then do it now via the magic 0xffffffff touch point with type 5 for sequence abort.
                                        if (!PsychHIDLastTouchEventFromEventBuffer(i, 0xffffffff)) {
                                            // Inject touch end event now ...
                                            PsychHIDAddEventToEventBuffer(i, &evt);

                                            // Then let the code below inject it again, but with type 5 for fail,
                                            // and magic id 0xffffffff:
                                            evt.type = 5;
                                            evt.rawEventCode = 0xffffffff;
                                        }
                                    }

                                    break;

                                case XI_TouchOwnership:
                                    // Ownership marker: We are the sole owner of this sequence, so got all the data
                                    // untampered :) - Set integrity bit for this touch point in last event for it:
                                    if (oldevt)
                                        oldevt->status |= (1 << 31);

                                    // printf("%i: XI_TouchOwnership!! %p\n", evt.rawEventCode, oldevt);

                                    break;
                            }

                            // Propagate positive integrity bit from preceeding event for this touch point:
                            if (oldevt && (oldevt->status & (1 << 31)))
                                evt.status |= (1 << 31);

                            // End of touch handling.
                        }

                        // Add anything but touch ownership events:
                        if (cookie->evtype != XI_TouchOwnership) {
                            // Update event buffer:
                            PsychHIDAddEventToEventBuffer(i, &evt);

                            // Tell waiting userspace (under KbQueueMutex protection for better scheduling) something interesting has changed:
                            PsychSignalCondition(&KbQueueCondition);
                        }
                    }

                    PsychUnlockMutex(&KbQueueMutex);
                }

                // Release event data:
                XFreeEventData(thread_dpy, cookie);
            }
        }

        // Done if we were only supposed to handle one event, which we did:
        if (blockingSinglepass) break;
    }

    return;
}

// Async processing thread for keyboard events:
void* KbQueueWorkerThreadMain(void* dummy)
{
    int rc;

    // Assign a name to ourselves, for debugging:
    PsychSetThreadName("PsychHIDKbQueue");

    // Try to raise our priority: We ask to switch ourselves (NULL) to priority class 2 aka
    // rt_fifo realtime scheduling, with a tweakPriority of +1, ie., raise the relative
    // priority level by +1 wrt. to the current level:
    if ((rc = PsychSetThreadPriority(NULL, 2, 1)) > 0) {
        printf("PsychHID: KbQueueStart: Failed to switch to realtime priority [%s].\n", strerror(rc));
    }

    while (1) {
        PsychLockMutex(&KbQueueMutex);

        // Check if we should terminate:
        if (KbQueueThreadTerminate) break;

        PsychUnlockMutex(&KbQueueMutex);

        // Perform X-Event processing until no more events are pending:
        KbQueueProcessEvents(TRUE);
    }

    // Done. Unlock the mutex:
    PsychUnlockMutex(&KbQueueMutex);

    // printf("DEBUG: THREAD TERMINATING...\n"); fflush(NULL);

    // Return and terminate:
    return(NULL);
}

int PsychHIDGetDefaultKbQueueDevice(void)
{
    int deviceIndex;
    XIDeviceInfo* dev = NULL;

    // Find first suitable slave device. For some reason, master keyboards don't work.

    // Whitelist scan: Use mouseemu virtual keyboard, if any:
    for(deviceIndex = 0; deviceIndex < ndevices; deviceIndex++) {
        dev = &info[deviceIndex];
        if ((dev->use == XISlaveKeyboard) && strstr(dev->name, "Mouseemu")) return(deviceIndex);
    }

    // Whitelist scan: Use first "*eyboard", if any. This is not foolproof, e.g., the
    // "Razer" gaming keyboard doesn't have the term "keyboard" in its name. But many
    // keyboards do:
    for(deviceIndex = 0; deviceIndex < ndevices; deviceIndex++) {
        dev = &info[deviceIndex];
        if ((dev->use == XISlaveKeyboard) && strstr(dev->name, "eyboard") && !strstr(dev->name, "XTEST")) return(deviceIndex);
    }

    // Blacklist scan: Use whatever comes first and isn't a button in disguise of
    // a keyboard or the virtual XTEST keyboard etc.:
    for(deviceIndex = 0; deviceIndex < ndevices; deviceIndex++) {
        dev = &info[deviceIndex];

        if ((dev->use == XISlaveKeyboard) && !strstr(dev->name, "XTEST") && !strstr(dev->name, "utton") && !strstr(dev->name, "Bus") &&
            !strstr(dev->name, "iSight") && !strstr(dev->name, "eceiver") && !strstr(dev->name, "amera") && !strstr(dev->name, "Gaming Mouse G502")) {
            return(deviceIndex);
        }
    }

    // Ok, no keyboard found at all. Fallback to mice, joysticks etc.:
    for(deviceIndex = 0; deviceIndex < ndevices; deviceIndex++) {
        dev = &info[deviceIndex];

        if ((dev->use == XISlavePointer) && !strstr(dev->name, "XTEST") && !strstr(dev->name, "utton") && !strstr(dev->name, "Bus") &&
            !strstr(dev->name, "iSight") && !strstr(dev->name, "eceiver") && !strstr(dev->name, "amera")) {
            return(deviceIndex);
        }
    }

    // Ok, nothing at all. Fallback to XTEST keyboard, which always exists. Maybe we deal with some
    // virtual keyboard like onscreen keyboard:
    for(deviceIndex = 0; deviceIndex < ndevices; deviceIndex++) {
        dev = &info[deviceIndex];
        if (dev->use == XISlaveKeyboard) return(deviceIndex);
    }

    // Nothing found? If so, abort:
    PsychErrorExitMsg(PsychError_user, "Could not find any useable keyboard device!");

    // Make compiler happy:
    return(-1);
}

PsychError PsychHIDOSKbQueueCreate(int deviceIndex, int numScankeys, int* scanKeys, int numValuators, int numSlots, unsigned int flags, unsigned int windowHandle)
{
    XIDeviceInfo* dev = NULL;

    // Valid number of keys?
    if (scanKeys && (numScankeys != 256)) {
        PsychErrorExitMsg(PsychError_user, "Second argument to KbQueueCreate must be a vector with 256 elements.");
    }

    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    } else if (deviceIndex >= ndevices) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. No such device!");
    }

    // Do we finally have a valid keyboard?
    dev = &info[deviceIndex];
    if (dev->use == XIMasterKeyboard) {
        PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. Master keyboards can not be handled by this function.");
    }

    // Keyboard queue for this deviceIndex already created?
    if (psychHIDKbQueueFirstPress[deviceIndex]) {
        // Yep. Release it, so we can start from scratch:
        PsychHIDOSKbQueueRelease(deviceIndex);
    }

    // Allocate and zero-init memory for tracking key presses and key releases:
    psychHIDKbQueueFirstPress[deviceIndex]   = calloc(256, sizeof(double));
    psychHIDKbQueueFirstRelease[deviceIndex] = calloc(256, sizeof(double));
    psychHIDKbQueueLastPress[deviceIndex]    = calloc(256, sizeof(double));
    psychHIDKbQueueLastRelease[deviceIndex]  = calloc(256, sizeof(double));
    psychHIDKbQueueScanKeys[deviceIndex]     = calloc(256, sizeof(int));
    memset(&psychHIDKbQueueOldEvent[deviceIndex], 0, sizeof(psychHIDKbQueueOldEvent[deviceIndex]));

    // Assign scanKeys vector, if any:
    if (scanKeys) {
        // Copy it:
        memcpy(psychHIDKbQueueScanKeys[deviceIndex], scanKeys, 256 * sizeof(int));
    } else {
        // None provided. Enable all keys by default:
        memset(psychHIDKbQueueScanKeys[deviceIndex], 1, 256 * sizeof(int));
    }

    // Store how many valuators we should return at most for this device:
    psychHIDKbQueueNumValuators[deviceIndex] = numValuators;

    // Store queue specific flags:
    psychHIDKbQueueFlags[deviceIndex] = flags;

    // Store associated X-Window handle, or zero for unspecified:
    psychHIDKbQueueXWindow[deviceIndex] = windowHandle;

    // Create event buffer:
    if (!PsychHIDCreateEventBuffer(deviceIndex, numValuators, numSlots)) {
        PsychHIDOSKbQueueRelease(deviceIndex);
        PsychErrorExitMsg(PsychError_system, "Failed to create keyboard queue due to out of memory condition.");
    }

    // Ready to use this keybord queue.
    return(PsychError_none);
}

void PsychHIDOSKbQueueRelease(int deviceIndex)
{
    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    }

    if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No such device!");
    }

    // Keyboard queue for this deviceIndex already exists?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Nothing to do then.
        return;
    }

    // Ok, we have a keyboard queue. Stop any operation on it first:
    PsychHIDOSKbQueueStop(deviceIndex);

    // Release its data structures:
    free(psychHIDKbQueueFirstPress[deviceIndex]); psychHIDKbQueueFirstPress[deviceIndex] = NULL;
    free(psychHIDKbQueueFirstRelease[deviceIndex]); psychHIDKbQueueFirstRelease[deviceIndex] = NULL;
    free(psychHIDKbQueueLastPress[deviceIndex]); psychHIDKbQueueLastPress[deviceIndex] = NULL;
    free(psychHIDKbQueueLastRelease[deviceIndex]); psychHIDKbQueueLastRelease[deviceIndex] = NULL;
    free(psychHIDKbQueueScanKeys[deviceIndex]); psychHIDKbQueueScanKeys[deviceIndex] = NULL;

    // Release kbqueue event buffer:
    PsychHIDDeleteEventBuffer(deviceIndex);

    // Done.
    return;
}

static void SingleXSelectEvents(XIEventMask *pemask, int deviceIndex, Window xwindow)
{
    Status rc;

    // Always grab a touch device for our exclusive use:
    if (XIMaskIsSet(pemask->mask, XI_TouchBegin) || XIMaskIsSet(pemask->mask, XI_RawTouchBegin)) {
        if (Success != (rc = XIGrabDevice(thread_dpy, pemask->deviceid, xwindow, CurrentTime, None, XIGrabModeAsync, XIGrabModeAsync, True, pemask)))
            printf("PsychHID-WARNING: KbQueueStart: Failed to grab touch input device %i with xinput device id %i for window %i: %s.\n",
                   deviceIndex, pemask->deviceid, (int) xwindow,
                   (rc == AlreadyGrabbed || rc == GrabFrozen) ? "Already grabbed by another application" :
                   (rc == GrabNotViewable) ? "Window not viewable" : "Unknown error");
    }
    else {
        // No touch device -> No exclusive grab:
        XISelectEvents(thread_dpy, xwindow, pemask, 1);
    }
}

// Helper: Set same event mask for all root windows of all x-screens for a server:
static void MultiXISelectEvents(XIEventMask *pemask, int deviceIndex, Window xwindow)
{
    int i;

    // xwindow zero ==> Apply to all X-Screens:
    if (xwindow == (Window) 0) {
        for (i = 0; i < ScreenCount(thread_dpy); i++)
            SingleXSelectEvents(pemask, deviceIndex, RootWindow(thread_dpy, i));
    }
    else if (xwindow <= (Window) ScreenCount(thread_dpy)) {
        // xwindow = [1 ; X-Screencount] ==> Apply to xwindow'th X-Screen:
        SingleXSelectEvents(pemask, deviceIndex, RootWindow(thread_dpy, ((int) xwindow) - 1));
    }
    else
        // Apply to specified xwindow handle:
        SingleXSelectEvents(pemask, deviceIndex, xwindow);
}

void PsychHIDOSKbQueueStop(int deviceIndex)
{
    psych_bool queueActive;
    int i;

    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    }

    if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No such device!");
    }

    // Keyboard queue for this deviceIndex already exists?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Nothing to do then.
        return;
    }

    // Keyboard queue already stopped?
    if (!psychHIDKbQueueActive[deviceIndex]) return;

    // Queue is active. Stop it:
    PsychLockMutex(&KbQueueMutex);

    // Setup event mask, so events from our associated xinput device
    // do not reach our event queue:
    XIEventMask emask;
    unsigned char mask[(XI_LASTEVENT + 7)/8];

    // Clear the event mask. That should to the trick:
    memset(mask, 0, sizeof(mask));
    emask.deviceid = info[deviceIndex].deviceid;
    emask.mask_len = sizeof(mask);
    emask.mask = mask;
    MultiXISelectEvents(&emask, deviceIndex, psychHIDKbQueueXWindow[deviceIndex]);
    XFlush(thread_dpy);

    // Mark queue logically stopped:
    psychHIDKbQueueActive[deviceIndex] = FALSE;

    PsychUnlockMutex(&KbQueueMutex);

    // Was this the last active queue?
    queueActive = FALSE;
    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
        queueActive |= psychHIDKbQueueActive[i];
    }

    // If more queues are active then we're done:
    if (queueActive) return;

    // No more active queues. Shutdown the common processing thread:
    PsychLockMutex(&KbQueueMutex);

    KbQueueThreadTerminate = TRUE;

    // Send some event to unblock the thread:
    XKeyEvent event;
    event.display = thread_dpy;
    event.window = DefaultRootWindow(thread_dpy);
    event.root = DefaultRootWindow(thread_dpy);
    event.subwindow = None;
    event.time = CurrentTime;
    event.x = 1;
    event.y = 1;
    event.x_root = 1;
    event.y_root = 1;
    event.same_screen = TRUE;
    event.type = KeyRelease;
    event.keycode = 0;
    event.state = 0;

    // printf("DEBUG: Sending termination fake event...\n"); fflush(NULL);
    XSelectInput(event.display, event.window, KeyReleaseMask);
    XFlush(thread_dpy);
    XSendEvent(event.display, event.window, TRUE, KeyReleaseMask, (XEvent *) &event);
    XFlush(thread_dpy);
    // printf("DEBUG: DONE.\n"); fflush(NULL);

    // Release grabbed touch input device:
    if ((psychHIDKbQueueNumValuators[deviceIndex] >= 4) && (PsychHIDIsTouchDevice(deviceIndex, NULL) >= 0))
        XIUngrabDevice(thread_dpy, info[deviceIndex].deviceid, CurrentTime);

    // Done.
    PsychUnlockMutex(&KbQueueMutex);

    // Shutdown the thread, wait for its termination:
    PsychDeleteThread(&KbQueueThread);
    KbQueueThreadTerminate = FALSE;

    // printf("DEBUG: THREAD JOINED.\n"); fflush(NULL);

    return;
}

void PsychHIDOSKbQueueStart(int deviceIndex)
{
    psych_bool queueActive;
    int i;
    int numValuators;

    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    }

    if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No such device!");
    }

    // Does Keyboard queue for this deviceIndex already exist?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Bad bad...
        printf("PsychHID-ERROR: Tried to start processing on non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No queue for that device yet!");
    }

    // Keyboard queue already stopped? Then we ain't nothing to do:
    if (psychHIDKbQueueActive[deviceIndex]) return;

    // Queue is inactive. Start it:

    // Will this be the first active queue, ie., aren't there any queues running so far?
    queueActive = FALSE;
    for (i = 0; i < PSYCH_HID_MAX_DEVICES; i++) {
        queueActive |= psychHIDKbQueueActive[i];
    }

    PsychLockMutex(&KbQueueMutex);

    // Clear out current state for this queue:
    memset(psychHIDKbQueueFirstPress[deviceIndex]   , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueFirstRelease[deviceIndex] , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueLastPress[deviceIndex]    , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueLastRelease[deviceIndex]  , 0, (256 * sizeof(double)));

    numValuators = psychHIDKbQueueNumValuators[deviceIndex];

    // Setup event mask, so events from our associated xinput device
    // get enqueued in our event queue:
    XIEventMask emask;
    unsigned char mask[(XI_LASTEVENT + 7)/8];

    memset(mask, 0, sizeof(mask));
    XISetMask(mask, XI_KeyPress);
    XISetMask(mask, XI_KeyRelease);

    // For mouse, joystick, gamepad and other "keyboardish" devices with buttons:
    XISetMask(mask, XI_RawButtonPress);
    XISetMask(mask, XI_RawButtonRelease);

    // Report motion events on pointing devices (mouse, trackpad etc.) if requested:
    // Queue flags +4 select for raw motion events, which don't go through pointer acceleration and other preprocessing:
    if ((numValuators >= 2) && (info[deviceIndex].use == XIMasterPointer || info[deviceIndex].use == XISlavePointer || info[deviceIndex].use == XIFloatingSlave))
        XISetMask(mask, (psychHIDKbQueueFlags[deviceIndex] & 0x4) ? XI_RawMotion : XI_Motion);

    // XInput version 2.2+ supported, and this device (Multi-)touch enabled?
    if ((numValuators >= 4) && (PsychHIDIsTouchDevice(deviceIndex, NULL) >= 0)) {
        // Due to what i suppose are X-Server bugs in touch-handling, normal touch
        // events don't work well at all on multi-X-Screen setups! The first (primary)
        // touch gets totally misscaled and wraps around -- some screwup in coordinate
        // transformations or state handling. This as of at least X-Server 1.19.6.
        //
        // Work around this: For single X-Screen setups, use regular touch events.
        // For multi X-Screen setups, use raw touch events instead, which behave more
        // reasonable, but have their own share of bugs -- easier managable bugs though...
        if (ScreenCount(thread_dpy) == 1) {
            // Single X-Screen: Regular touch events:
            XISetMask(mask, XI_TouchBegin);
            XISetMask(mask, XI_TouchUpdate);
            XISetMask(mask, XI_TouchEnd);
        }
        else {
            // Multi X-Screen: Raw events:
            XISetMask(mask, XI_RawTouchBegin);
            XISetMask(mask, XI_RawTouchUpdate);
            XISetMask(mask, XI_RawTouchEnd);
        }

        // Always request this, although it only seems to work for regular events,
        // and our touch sequence loss reporting is currently disabled due to its
        // own share of flaws...
        XISetMask(mask, XI_TouchOwnership);
    }

    emask.deviceid = info[deviceIndex].deviceid;
    emask.mask_len = sizeof(mask);
    emask.mask = mask;
    MultiXISelectEvents(&emask, deviceIndex, psychHIDKbQueueXWindow[deviceIndex]);
    XFlush(thread_dpy);

    // Mark this queue as logically started:
    psychHIDKbQueueActive[deviceIndex] = TRUE;

    // Queue started.
    PsychUnlockMutex(&KbQueueMutex);

    // If other queues are already active then we're done:
    if (queueActive) return;

    // No other active queues. We are the first one.

    // Start the common processing thread for all queues:
    PsychLockMutex(&KbQueueMutex);
    KbQueueThreadTerminate = FALSE;

    if (PsychCreateThread(&KbQueueThread, NULL, KbQueueWorkerThreadMain, NULL)) {
        // We are soo screwed:

        // Cleanup the mess:
        psychHIDKbQueueActive[deviceIndex] = FALSE;
        PsychUnlockMutex(&KbQueueMutex);

        // Whine a little bit:
        printf("PsychHID-ERROR: Start of keyboard queue processing failed!\n");
        PsychErrorExitMsg(PsychError_system, "Creation of keyboard queue background processing thread failed!");
    }

    // Up and running, we're done!
    PsychUnlockMutex(&KbQueueMutex);

    return;
}

void PsychHIDOSKbQueueFlush(int deviceIndex)
{
    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    }

    if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No such device!");
    }

    // Does Keyboard queue for this deviceIndex already exist?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Bad bad...
        printf("PsychHID-ERROR: Tried to flush non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No queue for that device yet!");
    }

    // Clear out current state for this queue:
    PsychLockMutex(&KbQueueMutex);
    memset(psychHIDKbQueueFirstPress[deviceIndex]   , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueFirstRelease[deviceIndex] , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueLastPress[deviceIndex]    , 0, (256 * sizeof(double)));
    memset(psychHIDKbQueueLastRelease[deviceIndex]  , 0, (256 * sizeof(double)));
    PsychUnlockMutex(&KbQueueMutex);

    return;
}

void PsychHIDOSKbQueueCheck(int deviceIndex)
{
    double *hasKeyBeenDownOutput, *firstPressTimeOutput, *firstReleaseTimeOutput, *lastPressTimeOutput, *lastReleaseTimeOutput;
    psych_bool isFirstPressSpecified, isFirstReleaseSpecified, isLastPressSpecified, isLastReleaseSpecified;
    int i;

    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    }

    if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No such device!");
    }

    // Does Keyboard queue for this deviceIndex already exist?
    if (NULL == psychHIDKbQueueFirstPress[deviceIndex]) {
        // No. Bad bad...
        printf("PsychHID-ERROR: Tried to check non-existent keyboard queue for deviceIndex %i! Call KbQueueCreate first!\n", deviceIndex);
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No queue for that device yet!");
    }

    // Allocate output
    PsychAllocOutDoubleArg(1, FALSE, &hasKeyBeenDownOutput);
    isFirstPressSpecified = PsychAllocOutDoubleMatArg(2, FALSE, 1, 256, 1, &firstPressTimeOutput);
    isFirstReleaseSpecified = PsychAllocOutDoubleMatArg(3, FALSE, 1, 256, 1, &firstReleaseTimeOutput);
    isLastPressSpecified = PsychAllocOutDoubleMatArg(4, FALSE, 1, 256, 1, &lastPressTimeOutput);
    isLastReleaseSpecified = PsychAllocOutDoubleMatArg(5, FALSE, 1, 256, 1, &lastReleaseTimeOutput);

    // Initialize output
    if(isFirstPressSpecified) memset((void*) firstPressTimeOutput, 0, sizeof(double) * 256);
    if(isFirstReleaseSpecified) memset((void*) firstReleaseTimeOutput, 0, sizeof(double) * 256);
    if(isLastPressSpecified) memset((void*) lastPressTimeOutput, 0, sizeof(double) * 256);
    if(isLastReleaseSpecified) memset((void*) lastReleaseTimeOutput, 0, sizeof(double) * 256);

    *hasKeyBeenDownOutput=0;

    // Compute and assign output:
    PsychLockMutex(&KbQueueMutex);

    for (i = 0; i < 256; i++) {
        double lastRelease  = psychHIDKbQueueLastRelease[deviceIndex][i];
        double lastPress    = psychHIDKbQueueLastPress[deviceIndex][i];
        double firstRelease = psychHIDKbQueueFirstRelease[deviceIndex][i];
        double firstPress   = psychHIDKbQueueFirstPress[deviceIndex][i];

        if (firstPress) {
            *hasKeyBeenDownOutput=1;
            if(isFirstPressSpecified) firstPressTimeOutput[i] = firstPress;
            psychHIDKbQueueFirstPress[deviceIndex][i] = 0;
        }

        if (firstRelease) {
            if(isFirstReleaseSpecified) firstReleaseTimeOutput[i] = firstRelease;
            psychHIDKbQueueFirstRelease[deviceIndex][i] = 0;
        }

        if (lastPress) {
            if(isLastPressSpecified) lastPressTimeOutput[i] = lastPress;
            psychHIDKbQueueLastPress[deviceIndex][i] = 0;
        }

        if (lastRelease) {
            if(isLastReleaseSpecified) lastReleaseTimeOutput[i] = lastRelease;
            psychHIDKbQueueLastRelease[deviceIndex][i] = 0;
        }
    }

    PsychUnlockMutex(&KbQueueMutex);

    return;
}

void PsychHIDOSKbTriggerWait(int deviceIndex, int numScankeys, int* scanKeys)
{
    int keyMask[256];
    int i;
    double t, tc;

    if (deviceIndex < 0) {
        deviceIndex = PsychHIDGetDefaultKbQueueDevice();
        // Ok, deviceIndex now contains our default keyboard to use - The first suitable keyboard.
    }

    if ((deviceIndex < 0) || (deviceIndex >= ndevices)) {
        // Out of range index:
        PsychErrorExitMsg(PsychError_user, "Invalid keyboard 'deviceIndex' specified. No such device!");
    }

    if(psychHIDKbQueueFirstPress[deviceIndex]) PsychErrorExitMsg(PsychError_user, "A queue for this device is already running, you must call KbQueueRelease() before invoking KbTriggerWait.");

    // Create a keyboard queue for this deviceIndex:
    memset(&keyMask[0], 0, sizeof(keyMask));
    for (i = 0; i < numScankeys; i++) {
        if (scanKeys[i] < 1 || scanKeys[i] > 256) PsychErrorExitMsg(PsychError_user, "Invalid entry for triggerKey specified. Not in valid range 1 - 256!");
        keyMask[scanKeys[i] - 1] = 1;
    }

    // Create keyboard queue with proper mask:
    PsychHIDOSKbQueueCreate(deviceIndex, 256, &keyMask[0], 0, 0, 0, 0);
    PsychHIDOSKbQueueStart(deviceIndex);

    PsychLockMutex(&KbQueueMutex);

    // Scan for trigger key:
    while (1) {
        // Wait until something changes in a keyboard queue:
        PsychWaitCondition(&KbQueueCondition, &KbQueueMutex);

        // Check if our queue had one of the dedicated trigger keys pressed:
        for (i = 0; i < numScankeys; i++) {
            // Break out of scan loop if key pressed:
            if (psychHIDKbQueueFirstPress[deviceIndex][scanKeys[i] - 1] != 0) break;
        }

        // Triggerkey pressed?
        if ((i < numScankeys) && (psychHIDKbQueueFirstPress[deviceIndex][scanKeys[i] - 1] != 0)) break;

        // No change for our trigger keys. Repeat scan loop.
    }

    // If we reach this point, we know some triggerkey has been pressed. As we aborted
    // the scan on detection of the first pressed key, we can't be certain we caught the
    // key with the earliest key press, maybe one of the untested keys was pressed even
    // earlier. Therefore do another pass over all keys to find the pressed one with the
    // earliest (minimum) pressed time:
    t = DBL_MAX;
    for (i = 0; i < numScankeys; i++) {
        tc = psychHIDKbQueueFirstPress[deviceIndex][scanKeys[i] - 1];
        if ((tc != 0) && (tc <= t)) t = tc;
    }

    // Done. Release the lock:
    PsychUnlockMutex(&KbQueueMutex);

    // Stop and release the queue:
    PsychHIDOSKbQueueStop(deviceIndex);
    PsychHIDOSKbQueueRelease(deviceIndex);

    // Return timestamp:
    PsychCopyOutDoubleArg(1, FALSE, t);

    return;
}
