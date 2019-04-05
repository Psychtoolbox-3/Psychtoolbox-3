# ppatest.py - Basic test of PsychPortAudio for Python
# Also tests WaitSecs, GetSecs, and PsychHID basics, as
# they are utilized here.
#
# Plays a beep tone, then another beep tone, then records sound for 10 seconds,
# then plays it back until the end, or until a key is pressed.
# The test will be subjec to change without notice, it is throwaway test code.
#
# (c) 2018 Mario Kleiner - Licensed under MIT license.

from psychtoolbox import *
from math import *
import numpy as np
#from soundfile import SoundFile

def printstatus(pahandle):
    info = PsychPortAudio('GetStatus', pahandle);
    print(info);
    return(info['Active']);

def run():

    # Samplerate:
    Fs = 48000;

    # Tone frequency:
    f = 500;

    # Length of audio vector - 1 secs of playback:
    sample = 1 * Fs;

    # Build sinewave:
    a = [0] * sample;

    for n in range(sample):
        a[n] = sin(2 * pi * f * n / Fs);

    # Replicated into two rows - One row for each stereo-channel, ie.
    # m'th row == m'th channel, n'th column = n'th sample frame:
    stereowav = [a, a];

    PsychPortAudio('Verbosity', 5);

    # Open audio device: Default audio device, default opmode (playback), default
    # latency/timing-precision (low-latency, high precision), sample rate Fs,
    # stereo (2) channel output:
    pahandle = PsychPortAudio('Open', [], [], [0], Fs, 2);
    PsychPortAudio('Volume', pahandle, 0.5);

    if True:
        # Fill in audio matrix for playback: columns = audio channels, rows = samples
        # Do it step-by-step for performance testing:
        stereowav = np.array(stereowav, order='f');
        # print('Type', type(stereowav), 'Shape', stereowav.shape, 'Datatype', stereowav.dtype, 'Order', stereowav.flags);
        stereowav = np.transpose(stereowav);
        # float32 is a tad faster than the default float64:
        stereowav = stereowav.astype('float32');
        # stereowav = np.zeros([10000,2], dtype='f')
    else:
        # Make it less boring:
        fname = '/home/kleinerm/Music/test.wav';
        myfile = SoundFile(fname);
        stereowav = myfile.read(dtype='float32', always_2d=True);
        #stereowav = myfile.read();
        myfile.close();

    print('Type', type(stereowav), 'Shape', stereowav.shape, 'Datatype', stereowav.dtype, 'Order', stereowav.flags);

    t1 = GetSecs();
    PsychPortAudio('FillBuffer', pahandle, stereowav);
    t2 = GetSecs();
    d1 = (1000 * (t2-t1));
    print('FillBuffer Duration', d1, ' msecs.');

    # Start playback for one repetition (1), 5 seconds from now, wait for sound onset:
    PsychPortAudio('Start', pahandle, 1, GetSecs() + 1, 1)

    # Go into a loop that prints playback status once a second, as long as playback
    # is active:
    info = PsychPortAudio('GetStatus', pahandle);
    print(info, ' Spec ', info['Active']);
    while printstatus(pahandle):
        WaitSecs('YieldSecs', 1);

    # Wait for sound to stop by itself, block until then:
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle, 1)
    print('StartTime', startTime, 'secs. Stop time', estStopTime, 'secs.\n');

    for n in range(sample):
        a[n] = sin(2 * pi * f * 1.5 * n / Fs);

    # Replicated into two rows - One row for each stereo-channel, ie.
    # m'th row == m'th channel, n'th column = n'th sample frame:
    stereowav2 = [a, a];
    stereowav2 = np.transpose(stereowav2);
    stereowav2 = stereowav2.astype('float32');

    t1 = GetSecs();
    b1 = PsychPortAudio('CreateBuffer', pahandle, stereowav2);
    t2 = GetSecs();
    d2 = (1000 * (t2-t1));
    print('CreateBuffer Duration', d2, ' msecs.');
    t1 = GetSecs();
    PsychPortAudio('FillBuffer', pahandle, b1);
    t2 = GetSecs();
    print('FillBuffer Duration', (1000 * (t2-t1)), ' msecs.');
    print('d2 / d1 = ', d2 / d1);

    PsychPortAudio('Start', pahandle, 1, GetSecs() + 1, 1)
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle, 1);
    print('StartTime', startTime, 'secs. Stop time', estStopTime, 'secs.\n');

    # Close sound device:
    PsychPortAudio('Close', pahandle);

    # 2nd test: Audio capture:
    pahandle = PsychPortAudio('Open', [], 2, 0, Fs, 2);

    # Preallocate an internal audio recording  buffer with a capacity of 10 seconds:
    PsychPortAudio('GetAudioData', pahandle, 10);

    # Start audio capture immediately and wait for the capture to start.
    # We set the number of 'repetitions' to zero,
    # i.e. record until recording is manually stopped.
    PsychPortAudio('Start', pahandle, 0, 0, 1);

    [audiodata, a, b, c] = PsychPortAudio('GetAudioData', pahandle, None, 9.9, None);
    print('Type', type(audiodata), 'Shape', audiodata.shape, 'Datatype', audiodata.dtype);

    PsychPortAudio('Close', pahandle);

    pahandle = PsychPortAudio('Open', [], [], [0], Fs, 2);

    PsychPortAudio('Volume', pahandle, 1);

    # Fill in audio matrix for playback:
    PsychPortAudio('FillBuffer', pahandle, audiodata);

    # Start playback for one repetition (1), 5 seconds from now, wait for sound onset:
    PsychPortAudio('Start', pahandle, 1, GetSecs() + 1, 1)

    # Go into a loop that prints playback status once a second, as long as playback
    # is active:
    info = PsychPortAudio('GetStatus', pahandle);
    print(info, ' Spec ', info['Active']);
    while printstatus(pahandle) and not PsychHID('KbCheck')[0]:
        WaitSecs('YieldSecs', 1);

    # Stop playback:
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle)
    print('StartTime', startTime, 'secs. Stop time', estStopTime, 'secs.\n');

    # Close sound device:
    PsychPortAudio('Close', pahandle);
