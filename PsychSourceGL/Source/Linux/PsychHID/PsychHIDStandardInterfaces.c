/*
  PsychToolbox3/Source/Linux/PsychHID/PsychHIDStandardInterfaces.c
  
  PROJECTS: PsychHID only.
  
  PLATFORMS:  Linux.
  
  AUTHORS:
  mario.kleiner@tuebingen.mpg.de    mk

  HISTORY:
  27.07.2011     mk     Created.
  
  TO DO:

*/

#include "PsychHIDStandardInterfaces.h"

static Display *dpy = NULL;
static Display *thread_dpy = NULL;
static int xi_opcode, event, error;
static XIDeviceInfo *info = NULL;
static int ndevices = 0;
static int masterDevice = -1;
static XDevice* x_dev[PSYCH_HID_MAX_DEVICES];

static	double* psychHIDKbQueueFirstPress[PSYCH_HID_MAX_DEVICES];
static	double* psychHIDKbQueueFirstRelease[PSYCH_HID_MAX_DEVICES];
static	double* psychHIDKbQueueLastPress[PSYCH_HID_MAX_DEVICES];
static	double* psychHIDKbQueueLastRelease[PSYCH_HID_MAX_DEVICES];
static  int*    psychHIDKbQueueScanKeys[PSYCH_HID_MAX_DEVICES];
static  psych_bool psychHIDKbQueueActive[PSYCH_HID_MAX_DEVICES];
static  psych_mutex KbQueueMutex;
static  psych_condition KbQueueCondition;
static	psych_bool  KbQueueThreadTerminate;
static  psych_thread KbQueueThread;
static	XEvent KbQueue_xevent;

static XDevice* GetXDevice(int deviceIndex)
{
	if (deviceIndex < 0 || deviceIndex >= PSYCH_HID_MAX_DEVICES) PsychErrorExitMsg(PsychError_user, "Invalid deviceIndex specified. No such device!");
	if (x_dev[deviceIndex] == NULL) x_dev[deviceIndex] = XOpenDevice(dpy, (XID) info[deviceIndex].deviceid);
	return(x_dev[deviceIndex]);
}

void PsychHIDInitializeHIDStandardInterfaces(void)
{
	int major, minor;
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

	// We must initialize XLib for multi-threaded operations / access on first
	// call:
	// TODO FIXME: We can only do this on Octave for now, not on Matlab!
	// Matlab uses XLib long before we get a chance to get here, but XInitThreads()
	// must be called as very first XLib function after process startup or bad things
	// will happen! So, we can't call it...
	// Because some system configurations can't handle multi-threaded x at all,
	// we allow users to opt-out of this if they define an environment variable
	// PSYCHTOOLBOX_SINGLETHREADEDX.
	#ifdef PTBOCTAVE3MEX
	if (NULL == getenv("PSYCHTOOLBOX_SINGLETHREADEDX")) XInitThreads();
	#endif

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

	// XInput V2 supported?
	major = 2;
	minor = 0;
	rc = XIQueryVersion(dpy, &major, &minor);
	if (rc == BadRequest) {
		printf("PsychHID: ERROR: No XI2 support. Server supports version %d.%d only.\n", major, minor);
		goto out;
	} else if (rc != Success) {
		printf("PsychHID: ERROR: Internal Error! This is a bug in Xlib.\n");
		goto out;
	}

	// printf("PsychHID: INFO: XI2 supported. Server provides version %d.%d.\n", major, minor);

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

	// XInput V2 supported?
	major = 2;
	minor = 0;
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

PsychError PsychHIDEnumerateHIDInputDevices(int deviceClass)
{
    const char *deviceFieldNames[]={"usagePageValue", "usageValue", "usageName", "index", "transport", "vendorID", "productID", "version", 
                                    "manufacturer", "product", "serialNumber", "locationID", "interfaceID", "totalElements", "features", "inputs", 
                                    "outputs", "collections", "axes", "buttons", "hats", "sliders", "dials", "wheels"};
    int numDeviceStructElements, numDeviceStructFieldNames=24, deviceIndex;
    PsychGenericScriptType	*deviceStruct;
    XIDeviceInfo *dev;
    int i, j;
    int numKeys, numAxis;
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
        PsychSetStructArrayDoubleElement("usagePageValue",	deviceIndex, 	(double) 1, deviceStruct);

	if (dev->use == XIMasterKeyboard || dev->use == XISlaveKeyboard) {
		// Usage 6 is for keyboard:
		PsychSetStructArrayDoubleElement("usageValue",	deviceIndex, (double) 6, deviceStruct);
	}
	else if (dev->use != XIFloatingSlave) {
		// Usage is 2 for mouse:
		PsychSetStructArrayDoubleElement("usageValue",	deviceIndex, (double) 2, deviceStruct);
	} else {
		// Assign meaningless Usage of zero for floating slave devices:
		PsychSetStructArrayDoubleElement("usageValue",	deviceIndex, (double) 0, deviceStruct);
	}

        PsychSetStructArrayStringElement("usageName",		deviceIndex, 	type, deviceStruct);
        PsychSetStructArrayDoubleElement("index",		deviceIndex, 	(double) i, deviceStruct);
        PsychSetStructArrayStringElement("transport",		deviceIndex, 	(dev->enabled) ? "enabled" : "disabled", deviceStruct);
        PsychSetStructArrayStringElement("product",		deviceIndex, 	dev->name, deviceStruct);
        PsychSetStructArrayDoubleElement("locationID",		deviceIndex, 	(double) dev->attachment, deviceStruct);
        PsychSetStructArrayDoubleElement("interfaceID",		deviceIndex, 	(double) dev->deviceid, deviceStruct);

        //PsychSetStructArrayDoubleElement("vendorID",		deviceIndex, 	(double)currentDevice->vendorID, 	deviceStruct);
        //PsychSetStructArrayDoubleElement("productID",		deviceIndex, 	(double)currentDevice->productID, 	deviceStruct);
        //PsychSetStructArrayDoubleElement("version",		deviceIndex, 	(double)currentDevice->version, 	deviceStruct);
        //PsychSetStructArrayStringElement("manufacturer",	deviceIndex, 	currentDevice->manufacturer, 		deviceStruct);
        //PsychSetStructArrayStringElement("serialNumber",	deviceIndex, 	currentDevice->serial, 			deviceStruct);

	numKeys = numAxis = 0;
	for (j = 0; j < dev->num_classes; j++) {
		if (dev->classes[j]->type == XIKeyClass) numKeys += (int) (((XIKeyClassInfo*) dev->classes[j])->num_keycodes);
		if (dev->classes[j]->type == XIButtonClass) numKeys += (int) (((XIButtonClassInfo*) dev->classes[j])->num_buttons);
		if (dev->classes[j]->type == XIValuatorClass) numAxis++;
	}

        PsychSetStructArrayDoubleElement("totalElements",	deviceIndex, 	(double) numKeys + numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("features",		deviceIndex, 	(double) dev->num_classes, deviceStruct);
        PsychSetStructArrayDoubleElement("inputs",		deviceIndex, 	(double) numKeys + numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("outputs",		deviceIndex, 	(double) 0, deviceStruct);
	PsychSetStructArrayDoubleElement("collections",         deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("axes",		deviceIndex, 	(double) numAxis, deviceStruct);
        PsychSetStructArrayDoubleElement("buttons",		deviceIndex, 	(double) numKeys, deviceStruct);
        PsychSetStructArrayDoubleElement("hats",		deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("sliders",		deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("dials",		deviceIndex, 	(double) 0, deviceStruct);
        PsychSetStructArrayDoubleElement("wheels",		deviceIndex, 	(double) 0, deviceStruct);
	
	deviceIndex++;
    }
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
		XQueryKeymap(dpy, keys_return);
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
		XQueryKeymap(dpy, keys_return);

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
	XIAnyClassInfo *classes;
	int i, j, dummy1, nclasses;

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
	double tnow;
	int i;
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

		if (FALSE) {
			if (KbQueue_xevent.type == KeyPress) {
				printf("KeyPress core event: key %i = %i\n", (int) ((XKeyEvent*) (&KbQueue_xevent))->keycode, (int) ((XKeyEvent*) (&KbQueue_xevent))->state);
				fflush(NULL);
			}
		}

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
				event = (XIDeviceEvent*) cookie->data;

				// Map Xinput device id to PTB 'deviceIndex' aka the proper keyboard queue:
				for (i = 0; i < ndevices; i++) if (event->deviceid == info[i].deviceid) break;

				// We're only interested in key press and release events, and only in
				// real ones, not XServer generated synthetic key auto-repeat events.
				// Also only for a device that we've registered:
				if ((i < ndevices) && !(event->flags & XIKeyRepeat) &&
				    ((cookie->evtype == XI_KeyPress) || (cookie->evtype == XI_KeyRelease) || (cookie->evtype == XI_ButtonPress) || (cookie->evtype == XI_ButtonRelease))) {

					// If this is a button event from a mouse/joystick etc. instead of a key event from a keyboard/keypad, then
					// all button indices are shifted by an offset of 1 for some weird reason. Decrement index by one to compensate
					// for this. [Tested on Ubuntu 10.10 and 11.10 with two mice and 1 joystick]
					if (((cookie->evtype == XI_ButtonPress) || (cookie->evtype == XI_ButtonRelease)) && (event->detail > 0)) event->detail--;

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

						key.keycode      = event->detail;
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
					if (psychHIDKbQueueActive[i] && (psychHIDKbQueueScanKeys[i][event->detail] != 0)) {
						// Yes: The queue wants to receive info about this key event.

						// Press or release?
						if ((cookie->evtype == XI_KeyPress) || (cookie->evtype == XI_ButtonPress)) {
							// Enqueue key press. Always in the "last press" array, because any
							// press at this time is the best candidate for the last press.
							// Only enqeue in "first press" if there wasn't any registered before,
							// ie., the slot is so far empty:
							if (psychHIDKbQueueFirstPress[i][event->detail] == 0) psychHIDKbQueueFirstPress[i][event->detail] = tnow;
							psychHIDKbQueueLastPress[i][event->detail] = tnow;
							evt.status |= (1 << 0);
						} else {
							// Enqueue key release. See logic above:
							if (psychHIDKbQueueFirstRelease[i][event->detail] == 0) psychHIDKbQueueFirstRelease[i][event->detail] = tnow;
							psychHIDKbQueueLastRelease[i][event->detail] = tnow;
							evt.status &= ~(1 << 0);
						}

						// Update event buffer:
						evt.timestamp = tnow;
						evt.rawEventCode = event->detail + 1;
						PsychHIDAddEventToEventBuffer(i, &evt);

						// Tell waiting userspace (under KbQueueMutex protection for better scheduling) something interesting has changed:
						PsychSignalCondition(&KbQueueCondition);
					}

					// Done with shared data access:
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

    // Blacklist scan: Use whatever comes first and isn't a button in disguise of
    // a keyboard or the virtual XTEST keyboard:
    for(deviceIndex = 0; deviceIndex < ndevices; deviceIndex++) {
        dev = &info[deviceIndex];
        if ((dev->use == XISlaveKeyboard) && !strstr(dev->name, "XTEST") && !strstr(dev->name, "Button") && !strstr(dev->name, "Bus")) return(deviceIndex);
    }
    
    // Nothing found? If so, abort:
    PsychErrorExitMsg(PsychError_user, "Could not find any useable keyboard device!");
}

PsychError PsychHIDOSKbQueueCreate(int deviceIndex, int numScankeys, int* scanKeys)
{
	XIDeviceInfo* dev = NULL;
	// int numKeys, j;

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
	if ((dev->use == XIMasterKeyboard) || (dev->use == XIMasterPointer)) {
		PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' specified. Master keyboards or master pointers cannot be handled by this function.");
	}

	/* This would detect the number of keys on the keyboard:
	numKeys = 0;
	for (j = 0; j < dev->num_classes; j++) {
		if (dev->classes[j]->type == XIKeyClass) numKeys += (int) (((XIKeyClassInfo*) dev->classes[j])->num_keycodes);
	}
	*/

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
    
	// Assign scanKeys vector, if any:
	if (scanKeys) {
		// Copy it:
		memcpy(psychHIDKbQueueScanKeys[deviceIndex], scanKeys, 256 * sizeof(int));
	} else {
		// None provided. Enable all keys by default:
		memset(psychHIDKbQueueScanKeys[deviceIndex], 1, 256 * sizeof(int));        
	}

	// Create event buffer:
	PsychHIDCreateEventBuffer(deviceIndex);

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

void PsychHIDOSKbQueueStop(int deviceIndex)
{
	psych_bool queueActive;
	int i;

	XIDeviceInfo* dev = NULL;

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
	XISelectEvents(thread_dpy, DefaultRootWindow(thread_dpy), &emask, 1);
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
	event.type = KeyPress;
	event.keycode = 0;
	event.state = 0;
	
	// printf("DEBUG: Sending termination fake event...\n"); fflush(NULL);
	XSelectInput(event.display, event.window, KeyPressMask);
	XFlush(thread_dpy);
	XSendEvent(event.display, event.window, TRUE, KeyPressMask, (XEvent *) &event);
	XFlush(thread_dpy);
	// printf("DEBUG: DONE.\n"); fflush(NULL);

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

	// Setup event mask, so events from our associated xinput device
	// get enqueued in our event queue:
	XIEventMask emask;
	unsigned char mask[(XI_LASTEVENT + 7)/8];

	memset(mask, 0, sizeof(mask));
	XISetMask(mask, XI_KeyPress);
	XISetMask(mask, XI_KeyRelease);

	// For mouse, joystick, gamepad and other "keyboardish" devices with buttons:
	XISetMask(mask, XI_ButtonPress);
	XISetMask(mask, XI_ButtonRelease);

	// XISetMask(mask, XI_Motion);

	emask.deviceid = info[deviceIndex].deviceid;
	emask.mask_len = sizeof(mask);
	emask.mask = mask;
	XISelectEvents(thread_dpy, DefaultRootWindow(thread_dpy), &emask, 1);
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
    PsychHIDOSKbQueueCreate(deviceIndex, 256, &keyMask[0]);
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
