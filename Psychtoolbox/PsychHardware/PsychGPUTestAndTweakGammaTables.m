function varargout = PsychGPUTestAndTweakGammaTables(win, xoffset, deviceType, injectFault, varargin)
%
% PsychGPUTestAndTweakGammaTables - Test and tweak GPU hardware gamma tables.
%
% This function is a helper function used for high-precision display
% devices like the Bits# from CRS and DataPixx/ViewPixx/ProPixx from VPixx,
% which need identity passthrough of framebuffer pixel data to the video outputs
% without any interference by the GPU and other intermediate encoder circuitry.
%
% It augments functions like LoadIdentityClut, trying to cope with GPUs that
% are too broken to be fixable just by that function without any actual
% measurement of video signals.
% The function is called by toolboxes for special display output hardware
% and makes use of / calls back into hardware specific functions of that
% hardware to allow an iterative feedback loop of parameter tweaks and
% measurements to execute in order to optimize gamma tables for optimal
% pixel passthrough after other hardware-independent measures have failed.
%
% The routine is, e.g., used by the 'GPUEncoderTest' / 'CheckGPUSanity'
% functions of the PsychDatapixx and BitsPlusPlus functions.
%
% Input parameters: (Mandatory)
% 'win' = Onscreen window handle for onscreen window which displays to the
%         external high precision display device.
%
% 'deviceType' = Type of device: 0 = VPixx Inc. Data-/View-/Pro-Pixx.
%                                1 = CRS Bits#
% 
% 'injectFault' = 1 = Intentionally setup a slightly faulty LUT to perturb
%                 the signal and test the tweaking procedure. 0 = Don't.
%
% Returns 0 on success, 1 on failure.
%

% History:
% 11/04/2013  mk  Written. Consolidated from PsychDataPixx driver function
%                 to share all the common code for VPixx and CRS products.
%

global GL;

    if nargin < 4
        error('Required minimum 3 parameters missing!');
    end
    
    if deviceType == 0
        % Open our own connection to Datapixx:
        devName = 'PsychDataPixx';
        PsychDataPixx('Open');
        
        % Repeat each measurement 60 times on DataPixx:
        numRescans = 60;
    end

    if deviceType == 1
        devName = 'BitsPlusPlus';

        % Retrieve IOPort handle to Bits# serial port connection:
        bitsSharpPort = varargin{1};

        % Repeat each measurement 20 times on Bits# to compensate for its
        % slow readback:
        numRescans = 20;
        
        % Switch to Bits# status screen before pixel readback. Why? Pixel readback would
        % implicitely switch to status screen if it is not already active, but this
        % switch takes multiple seconds and causes instable operation of pixel readback,
        % ie., it delays readback and can cause timeouts in our readback code.
        % We prevent this by switching beforehand and giving it more than enough time to
        % stabilize:
        fprintf('%s: INFO: Switching Bits# into diagnostic mode. Will take about 5 seconds...\n', devName);
        IOPort('Write', bitsSharpPort, ['$statusScreen' char(13)]);
        WaitSecs('YieldSecs', 5);
    end

    % Set to 1 for fault injection to test the procedure:
    if injectFault
        fprintf('%s: INFO: Injecting perturbed gamma table to self-test my tweaking abilities...\n', devName);
        oldlut = Screen('ReadNormalizedGammaTable', win);
        Screen('LoadNormalizedGammaTable', win, oldlut * 0.95);
        %Screen('LoadNormalizedGammaTable', win, rand(size(oldlut)));
    end
    
    % Build 256 pixel level RGB test gradient:
    testdata = uint8(zeros(3, 256));
    testdata(:,:) = uint8(repmat(0:255, 3, 1));

    % Clear the onscreen window to black and set its clear color to black.
    Screen('FillRect', win, 0);

    % Wait at least two video refresh cycles to make sure the LUT's in the
    % GPU have stabilized:
    if deviceType == 0
        % DataPixx: Wait for two vsync pulses:
        Datapixx('RegWrRdVideoSync');
        Datapixx('RegWrRdVideoSync');
    else
        % Other: Just wait for approximately 2 refresh durations:
        WaitSecs('YieldSecs', 2 * Screen('GetFlipInterval', win));
    end

    % At this point the GPU gamma tables should be already loaded with an
    % identity LUT via LoadIdentityClut(). The GL struct and moglcore are
    % initialized so we can use low-level glXXX commands.
    
    % Testing loop: Will display the horizontal test gradient in the top
    % scanline of the display, then read it back from the device. If
    % everything works correctly then the readback gradient should match
    % the drawn gradient. A mismatch means that either the gamma tables are
    % not loaded with a proper identity LUT, or that some kind of
    % spatio-temporal display dithering is active and interfering.
    %
    % We can try to auto-correct for wrong identity LUT's with our iterative
    % correction loop. It will tweak the gamma tables until a match is
    % achieved. We can't do anything about dithering though beyond
    % reporting it:
    oldlut = Screen('ReadNormalizedGammaTable', win);
    nrlutslots = size(oldlut, 1);
    curlut = oldlut;
    
    psychlasterror('reset');

    % We tweak gamma table elements in steps of 1/1024 units per iteration.
    % As the granularity of Digital DVI-D is only 8 bits or 1/256 units,
    % that means we are resolving 4 times more fine-grained than the device
    % resolution on DVI. This just for safety to avoid "jumping into the
    % wrong bin" due to roundoff errors in the OS/Driver. Could probably do
    % more coarse-grained, but better safe than sorry...
    stepsize = 1/1024;
    nIter = 100;
    failcount = 0;
    successcount = 0;
    regressioncount = 0;
    fluctcount = 0;
    retry = 1;
    
    % Draw test pattern:
    glRasterPos2i(xoffset, 10);
    glDrawPixels(256, 10, GL.RGB, GL.UNSIGNED_BYTE, repmat(testdata, 1, 10));

    % Show it at next retrace:
    Screen('Flip', win);
    
    % Wait generous 5 retrace cycles in case a compositor is delaying things:
    Screen('WaitBlanking', win, 5);
    
    fprintf('%s: INFO: GPU gamma table and display encoder test running. Please be patient...\n', devName);
    tStart = GetSecs;
    
    % Enter test and auto-correction loop:
    while retry
        retry2=1;
        while retry2 && (fluctcount < nIter)
            if deviceType == 0
                % Wait for next vsync pulse:
                Datapixx('RegWrRdVideoSync');

                % Read back stabilized pixel data:
                oldrealpixels = Datapixx('GetVideoLine', 256);
            end

            if deviceType == 1
                % Wait for next vsync:
                Screen('WaitBlanking', win);

                % Readback 2nd topmost scanline from Bits# :
                oldrealpixels = BitsPlusPlus('GetVideoLine', 256, 2);
            end

            % Repeat measurement until at least numRescans consecutive scans give
            % the same value:
            for i=1:numRescans
                if deviceType == 0
                    % Wait for next vsync pulse:
                    Datapixx('RegWrRdVideoSync');

                    % Read back stabilized pixel data:
                    realpixels = Datapixx('GetVideoLine', 256);
                end

                if deviceType == 1
                    % Wait for next vsync:
                    Screen('WaitBlanking', win);

                    % Readback 2nd topmost scanline from Bits# :
                    realpixels = BitsPlusPlus('GetVideoLine', 256, 2);
                end

                if any(any(int16(oldrealpixels) - int16(realpixels)))
                    % Difference between 1st measurement and i'th
                    % measurement detected. Abort scan and retry:
                    fluctcount = fluctcount + 1;
                    clc;
                    fprintf('%s: INFO: GPU gamma table and display encoder test running. Please be patient...\n', devName);
                    fprintf('Instable DVI-D signal, possibly due to dithering! Will try to fix this by iterative tweaking of gamma tables.\n');
                    fprintf('Iteration Nr. %i in progress. Current pixelvalues received from topmost scanline:\n\n', failcount);
                    disp(realpixels);
                    break;
                else
                    if i==numRescans
                        retry2=0;
                    end
                end
            end
        end
        
        if fluctcount >= nIter
            % Totally unstable DVI signal. This could be due to
            % (spatio-)temporal dithering being enabled on the display head
            % or due to serious signal integrity issues.
            fprintf('%s: INFO: No stable signal received during over %i iterations! Interference by dithering?!?\n', devName, nIter);
            fluctcount = 0;     
        end
        
        % 'realpixels' contains a measurement that has been stable for at
        % least numRescans video cycles. We use this to judge the quality of the
        % current GPU LUT settings and perform corrections if needed.
        
        % Compute delta matrix to reference values: This should contain
        % all-zeros on properly working GPU's:
        deltavec = int16(testdata) - int16(realpixels);
        
        % Difference?
        if any(any(deltavec))
            % Ohoh!
            if failcount == 0
                fprintf('%s: WARNING: Your GPU has a wrong identity gamma table loaded, possibly due to driver bugs.\n', devName);
                fprintf('%s: WARNING: I will try to automatically compensate for this now. Please standby...\n', devName);
            end
            
            failcount = failcount + 1;
            
            if successcount > 0
                regressioncount = regressioncount + 1;
                fprintf('%s: WARNING: %i. regression detected after tweak %i! Previously successfull setting produces wrong values now\nafter %i successfull iterations!\n', devName, regressioncount, failcount, successcount);
                successcount = 0;
            end
            
            % Convert to double matrix and transpose to match LUT format.
            % Then convert to a corrective increment of stepsize units per
            % 1 pixel unit error:
            deltavec = transpose(double(deltavec)) * stepsize;
            
            % Build correction vector tweakvec: deltavec is a 256 rows by 3
            % columns matrix of delta values to add to curlut for tweaking
            % it to a more correct value. curlut, as read from the GPU,
            % always has 3 columns for (Red, Green, Blue) just like
            % deltavec, but it can have more than 256 slots on high-end
            % GPU's with bigger gamma tables. E.g., the NVidia
            % QuadroFX-3800 when run under Linux with the binary blob
            % driver will have a 2^11 = 2048 slot lut. Therefore we need to
            % "stretch" deltavec vertically so it matches the size of the
            % 'curlut'. By design of the Datapixx/Bits# and DVI-D data
            % transmission, we can always only compute 256 rows of
            % meaningful tweaking information, therefore we simply repmat
            % replicate successive rows in deltavec to blow it up to a
            % matching tweakvec. This should be good enough to init the 256
            % rows subset of curlut we actually care about with the tweaked
            % values we found:
            tweakvec = zeros(size(curlut));
            chunksize = nrlutslots / 256;
            for tchunk = 0:255
                tweakvec((1 + (tchunk * chunksize)):((tchunk * chunksize) + chunksize), :) = repmat(deltavec(1 + tchunk, :), chunksize, 1);
            end
            
            % Update gamma lut:
            curlut = curlut + tweakvec;
            
            mmcount = sum(sum(deltavec ~= 0));
            fprintf('%s: INFO: Gamma table tweak iteration %i. Had %i mismatching pixels. Retesting with updated LUT...\n', devName, failcount, mmcount);
            
            % Clamp to valid 0-1 range:
            curlut(curlut > 1) = 1;
            curlut(curlut < 0) = 0;
            
            % Upload modified LUT to GPU: We try-catch protect this to
            % avoid error-abort on MS-Windows if invalid gamma tables are
            % encountered:
            try
                Screen('LoadNormalizedGammaTable', win, curlut);
            catch
                psychlasterror('reset');
            end
            
            % Tried too long, too hard?
            if failcount > 255
                % Running with non-zero rasterizer offset?
                if xoffset ~= 0
                    % Hmm. Could be a cascade of bugs. Reset xoffset to
                    % zero and retry whole procedure just to be on the safe
                    % side:
                    fprintf('%s: WARNING: No success so far with rasterizeroffset %i. Retrying with zero-offset...\n', devName, xoffset);
                    failcount = 0;
                    curlut = oldlut;

                    Screen('LoadNormalizedGammaTable', win, oldlut);
                    
                    % Draw test pattern:
                    Screen('Flip', win);
                    xoffset = 0;
                    glRasterPos2i(xoffset, 10);
                    glDrawPixels(256, 10, GL.RGB, GL.UNSIGNED_BYTE, repmat(testdata, 1, 10));

                    % Show it at next retrace:
                    Screen('Flip', win);
                    Screen('WaitBlanking', win);
  
                    % Encore une fois...
                    continue;
                end
                
                % Ok, this can't be. Restore gamma table and give up:
                Screen('LoadNormalizedGammaTable', win, oldlut);
                                
                % Fail:
                fprintf('%s: ERROR: Did not reach a stable result after 256 tweak iterations.\n', devName);
                fprintf('%s: ERROR: This is hopeless, something is seriously wrong with your GPU. I give up!\n', devName);
                fprintf('%s: ERROR: One possible reason could be that (spatio-)temporal dithering is enabled for your\n', devName);
                fprintf('%s: ERROR: display due to some graphics driver bug, which prevents me from converging\n', devName);
                fprintf('%s: ERROR: to a stable result. Your high precision device will not be useable for high-precision\n', devName);
                fprintf('%s: ERROR: color or luminance displays with this broken graphics card/driver. Fix your system!\n', devName);
                
                if deviceType == 0
                    % Close our own connection to Datapixx:
                    PsychDataPixx('Close');
                end

                varargout{1} = 1;
                return;
            end
        else
            % No difference. Good. Let's repeat the exercise a few times to
            % rule out bad interference from temporal dithering:
            successcount = successcount + 1;

            fprintf('%s: INFO: Posttest iteration %i passed. Good! [Total elapsed time so far %f seconds.]\n', devName, successcount, GetSecs - tStart);
            
            if successcount >= 10
                % 10 successfull consecutive repetitions. That's good
                % enough for us:
                retry = 0;
            end
        end
    end

    % Ok, we managed to get through:
    Screen('Flip', win);
    
    % How did we do?
    if (successcount >= 10) && (failcount == 0)
        % Purrfect! Return full success:
        fprintf('%s: INFO: Online test of GPU hardware identity clut and display encoder passed after a total of %f seconds.\n', devName, GetSecs - tStart);
        varargout{1} = 0;
    else
        % Got through, but with a non-zero failcount. That probably means
        % that display temporal dithering is not an issue, but that the
        % gamma table needed tweaking.
        fprintf('%s: INFO: Managed to find a useable gamma table setting after %i tweak iterations and %f seconds. :-)\n', devName, failcount, GetSecs - tStart);
        if regressioncount > 0
            fprintf('%s: INFO: There were %i regressions throughout the procedure. Could be problems with active display dithering or bad signal quality?\n', devName, regressioncount);
        end
        if fluctcount == 0
            fprintf('\n%s: INFO: There wasn''t any instability during the last testing cycle,\n', devName);
            fprintf('%s: INFO: so harmful display dithering may be absent now and you will be probably fine.\n', devName);
        end
        fprintf('%s: INFO: I will store this LUT on your filesystem for future use.\n', devName);
        if 0
            fprintf('%s: INFO: The difference matrix between the original LUT and the corrected LUT looks like this:\n\n', devName);
            disp(curlut - oldlut);
        end
        fprintf('\n\n');
        
        % Store the current gamma table LUT in a config file for later use
        % by LoadIdentityClut in successive runs:
        SaveIdentityClut(win);
        varargout{1} = 0;
    end

    if deviceType == 0
        % Close our own connection to Datapixx:
        PsychDataPixx('Close');
    end

    return;
end
