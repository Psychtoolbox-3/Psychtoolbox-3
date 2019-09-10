"""

PsychHID('KbQueueCreate', [deviceNumber][, keyFlags=all][, numValuators=0][, numSlots=10000][, flags=0][, windowHandle=0])
PsychHID('KbQueueRelease' [, deviceIndex])
[navail] = PsychHID('KbQueueFlush' [, deviceIndex][, flushType=1])
PsychHID('KbQueueStart' [, deviceIndex])
PsychHID('KbQueueStop' [, deviceIndex])
[keyIsDown, firstKeyPressTimes, firstKeyReleaseTimes, lastKeyPressTimes, lastKeyReleaseTimes]=PsychHID('KbQueueCheck' [, deviceIndex])
secs=PsychHID('KbTriggerWait', KeysUsage, [deviceNumber])
[event, navail] = PsychHID('KbQueueGetEvent' [, deviceIndex][, maxWaitTimeSecs=0])


Support for access to generic USB devices: See 'help ColorCal2' for one usage example:


usbHandle = PsychHID('OpenUSBDevice', vendorID, deviceID [, configurationId=0])
PsychHID('CloseUSBDevice' [, usbHandle])
outData = PsychHID('USBControlTransfer', usbHandle, bmRequestType, bRequest, wValue, wIndex, wLength, inData)

"""
import sys
from . import PsychHID


def num_devices():
    """Returns the number of connected HID devices"""
    return PsychHID('NumDevices')


def devices(device_class=None):
    """Returns a list of devices of a given class

    :param device_class: integer
    :return: list of devices
    """
    return PsychHID('Devices', device_class)


def get_keyboard_indices(name='', serial_number=''):
    """This is the equivalent of PTB GetKeyboardIndices function

    :param name: string (optional)
        name to match
    :param serial_number: string  (optional)
        serial number to match
    :return: (keyboardIndices, productNames, allInfos)
    """
    keyboardIndices = []
    productNames = []
    allInfos = []
    if sys.platform != 'darwin':
        devs = devices(device_class=4)  # on win/linux this is just keyboards
    else:
        devs = devices()
    for ii, dev in enumerate(devs):
        # filter out non-matches
        # print("{}: {}, {}, {}".format(
        #     dev['product'], dev['usageValue'], dev['usagePageValue'], dev['serialNumber']))
        if dev['usagePageValue'] != 1 or dev['usageValue'] != 6:
            continue  # wrong spec - doesn't look like a keyboard
        if name and name not in dev['product']:
            continue
        if serial_number and serial_number != dev['serialNumber']:
            continue  # doesn't match the serial number
        # we've got a match so store the data
        # for mac we could also check if 'Keyboard' in dev['usageName'] but not win32
        keyboardIndices.append(dev['index'])
        productNames.append(dev['product'])
        allInfos.append(dev)
    return keyboardIndices, productNames, allInfos


class Device:
    def __init__(self, device_number):
        self.device_number = device_number

    def num_elements(self):
        return PsychHID('NumElements', self.device_number)

    def elements(self):
        return PsychHID('Elements', self.device_number)

    def num_collections(self):
        return PsychHID('NumCollections', self.device_number)

    def collections(self):
        return PsychHID('Collections', self.device_number)

    def get_state(self, element_number, calibrated=False):
        if calibrated:
            return PsychHID('CalibratedState', self.device_number,
                            element_number)
        else:
            return PsychHID('RawState', self.device_number, element_number)

    def get_report(self, report_type, report_id, report_bytes):
        """Returns a report and error code for the device

        :param report_type:
        :param report_id:
        :param report_bytes:
        :return: [report, err]
        """
        return PsychHID('GetReport', self.device_number)

    def set_report(self, report_type, report_id, report):
        """"""
        return PsychHID('SetReport', self.device_number,
                        report_type, report_id, report)

    def give_me_reports(self, report_bytes=None):
        """:returns: [reports, err]"""
        return PsychHID('GiveMeReports', self.device_number, report_bytes)

    def receive_reports(self, options=None):
        """:returns: err code"""
        return PsychHID('ReceiveReports', self.device_number, options)

    def stop_reports(self):
        """:returns: err code"""
        return PsychHID('ReceiveReportsStop', self.device_number)


class Keyboard():
    def __init__(self, device_number=None, buffer_size=10000):
        """A Keyboard object is like a Device() with key-specific functions

        :param device_number: float or int

            Can be found from

        :param buffer_size: float or int

            Determines how many key events (up and down) are stored
        """
        self.device_number = device_number
        self._create_queue(buffer_size)

    def check(self, scan_list=None):
        """Checks for events """
        return PsychHID('KbCheck', self.device_id, scan_list)

    def _create_queue(self, num_slots=10000, flags=0, win_handle=0):
        PsychHID('KbQueueCreate', self.device_number,
                 None, 0, num_slots, flags, win_handle)
        # [deviceNumber][, keyFlags=all][, numValuators=0][, numSlots=10000]
        # [, flags=0][, windowHandle=0])

    def _release_queue(self):
        PsychHID('KbQueueRelease', self.device_number)

    def flush(self, flush_type=1):
        """Flushes the keybard queue and returns the number of available evts"""
        return PsychHID('KbQueueFlush', self.device_number, flush_type)

    def queue_start(self):
        """Starts recording key presses to the queue"""
        PsychHID('KbQueueStart', self.device_number)

    def queue_stop(self):
        """Stops recording key presses to the queue"""
        PsychHID('KbQueueStop', self.device_number)

    def queue_check(self):
        """Returns key events from the buffer
        :param self:
        :return: (keyIsDown,
            firstKeyPressTimes, firstKeyReleaseTimes,
            lastKeyPressTimes, lastKeyReleaseTimes)
        """
        return PsychHID('KbQueueCheck', self.device_number)

    def queue_get_event(self, max_wait_secs=0):
        """

        :param max_wait_secs: float
        :return: (event, navail)
        """
        return PsychHID('KbQueueGetEvent', self.device_number, max_wait_secs)

    def trigger_wait(self, keys):
        secs = PsychHID('KbTriggerWait', keys, self.device_number)

    def start_trapping(self):
        PsychHID('Keyboardhelper', -12)

    def stop_trapping(self):
        PsychHID('Keyboardhelper', -10)
