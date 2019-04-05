# hidtest.py - Basic test of PsychHID for Python
# Also tests WaitSecs, as it is utilized here.
#
# Collects key responses for 5 seconds via KbQueues while being blocked on
# the main thread. The disabled latter part will perform some low-level access
# to a USB-HID device, to test HID report send/receive. However, this is so
# device specific as being useless without customization to your specific machine
# OS and HID device. Rest assured it has some use for Mario's laptop.
#
# The test will be subjec to change without notice, it is throwaway wip code.
#
# (c) 2018 Mario Kleiner - Licensed under MIT license.

from psychtoolbox import *

def run():

    # KbQueueTest
    WaitSecs('YieldSecs', 1);
    keys = [1] * 256;
    PsychHID('KbQueueCreate', [], keys);
    PsychHID('KbQueueStart');
    PsychHID('Keyboardhelper', -12);
    WaitSecs('YieldSecs', 5);
    PsychHID('KbQueueStop');
    PsychHID('Keyboardhelper', -10);
    while PsychHID('KbQueueFlush'):
        evt = PsychHID('KbQueueGetEvent');
        print(evt);

    PsychHID('KbQueueRelease');
    return;

    # USB HID reports test:
    #hiddevs = PsychHID('Devices')
    mouseIndex = 5;
    wheelDelta = 0;
    WaitSecs('YieldSecs', 1);

    options = dict();
    options['print'] = 0;
    err = PsychHID('ReceiveReports', mouseIndex, options);
    print(err);

    while not PsychHID('KbCheck')[0]:
        rep = PsychHID('GetReport', mouseIndex, 1, 0, 5);
        print(rep);
        WaitSecs('YieldSecs', 1);

    err = PsychHID('ReceiveReportsStop', mouseIndex)
    print(err);

    print('\n\nNow all the other remaining ones:\n\n');
    [rep, errs] = PsychHID('GiveMeReports', mouseIndex);
    print(errs);
    print('\n\nReports:\n\n');
    print(rep);
