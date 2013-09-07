function HDRViewer(imfilepattern, dummymode, sf)
% HDRViewer([imfilepattern][, dummymode][, scalefactor]) -- Load and show high
% dynamic range images on the BrightSide Technologies High Dynamic Range
% display device.
%
% The viewer allows to cycle through a sequence of images in the given
% folder, matching the given filename pattern. It allows to adjust the
% intensity online, as well as to zoom into regions of the image.
%
% 'imfilepattern' - Filename search pattern of the HDR images to load,
% e.g., 'myimages*.hdr' would load all images starting with 'myimages' and
% ending with extension '.hdr'. 
%
% 'dummymode' - If set to 1 we only run in emulation mode without use of
% the HDR device or BrightSide core library.
%
% Control keys:
%
% You can cycle through all images by pressing the space-bar. You can change
% image intensity scaling by pressing the cursor up- and down-keys. You can
% exit the viewer by pressing ESCape. By pressing the 'q' key, you can
% toggle between float precision HDR and 256 level 8 bit quantized output.
%
% Zoom mode:
% Pressing the left mouse button will enable zoom mode. Dragging the mouse
% while the button is depressed allows to change the region which should be
% zoomed. Releasing the mouse button will show the selected region zoomed
% up to fullscreen. A single mouse click will reset to the full image.
%
% Btw. you can find lot's of nice free HDR images by Googling on the
% internet!
%
% This is an extension of SimpleHDRDemo.m, see that file for basic usage of
% HDR displays with Psychtoolbox.
%
% History:
% Written 2006 by Mario Kleiner - MPI for Biological Cybernetics, Tuebingen, Germany
% and Oguz Ahmet Akyuz - Department of Computer Science, University of Central Florida.

% Make sure we run on OpenGL-Psychtoolbox. Abort otherwise.
PsychDefaultSetup(1);

% Shorten display timing tests - Pretty futile on such a slow device
% anyway...
Screen('Preference', 'SkipSyncTests', 1);

% Startup with black screen to not hurt our eyes:
Screen('Preference', 'VisualDebugLevel', 3);

% Run demo in dummy mode? We default to the real thing.
if nargin < 2 || isempty(dummymode)
    dummymode = 0;
end

% Name of an image file passed? If not, we default to all *.hdr files in
% current working directory:
if nargin < 1 || isempty(imfilepattern)
    % No pattern provided. Search in current working directory.
    imfilepattern = '*.hdr';
end

if nargin < 3
    sf = [];
end

imfilenames = dir(imfilepattern)

try
    % Disable keypress output to Matlab window:
    ListenChar(2);
    
    % Find screen to display: We choose the one with the highest number,
    % assuming this is the HDR display:
    screenid=max(Screen('Screens'));
    
    % Open a standard fullscreen onscreen window, double-buffered with black
    % background color instead of the default white one: win is the window
    % handle for this window. We use the BrightSideHDR() command instead of
    % Screen(). It is a convenience wrapper around Screen, doing all the
    % additional setup work for the HDR display:
    if dummymode
        % Dummy mode: Don't run on real HDR display:
        win = BrightSideHDR('DummyOpenWindow', screenid, 0);
    else
        % HDR mode: Setup everything for HDR rendering:
        win = BrightSideHDR('OpenWindow', screenid, 0);
    end

    winrect = Screen('Rect', win);
    
    iteration = 0;
    abortit=0;

    % Set text properties:
    Screen('TextSize', win, 24);
    Screen('TextStyle', win , 0);

    % Multiply all image values by a scaling factor: The HDR accepts
    % values between zero and infinity, but the useable range seems to be
    % zero (Dark) to something around 3000-4000. At higher values, it
    % saturates in a non-linear fashion. Fractional values, e.g,. 0.5 are
    % resolved at an unknown quantization level, so the range of
    % displayable intensity levels is more than 3000-4000 steps. We don't
    % know the real resolution without proper calibration, but according to
    % their Siggraph 2004 paper its supposed to be more than 14000 levels,
    % covering the full operating range in steps of single JND's. Who
    % knows...

    % Initial flip to black background:
    Screen('Flip', win);
    
    while ~abortit
        iteration=iteration + 1;
        imagename=imfilenames(mod(iteration, size(imfilenames,1))+1).name;
        
        % Show some status info:
        DrawFormattedText(win, ['Loading image: ' imagename '...'], 'center', 'center', [0 255 0]);
        Screen('Flip', win);
        
        % Read image into Matlab matrix:
        img = HDRRead(imagename, 1);
        if isempty(img)
            continue;
        end

        % Set step size for change of scaling factor:
        sfstep = 3500 / max(max(max(img))) * 0.01;
        
        % Scale intensities by some factor. Here we set the default value:
        if isempty(sf)
            if ~dummymode
                sf = 3500 / max(max(max(img)));
            else
                sf = 0.005;
            end
        end
        
        % Build a Psychtoolbox 32 bpc float texture from the image array
        % by setting the (optional) 'floatprecision' flag to 2.
        texid = Screen('MakeTexture', win, img, [], 2, 2);

        % Build also a version of the image that is quantized to 8 Bit:
        quantimg = uint8((img / max(max(max(img)))) * 255);
        ldrtexid = Screen('MakeTexture', win, quantimg, [], 2);
        
        needupdate = 1;
        zoomset = 0;
        xleft=0;
        xtop=0;
        xright=-1;
        xbottom=-1;
        hdrmode = 1;
        
        while(needupdate)            
            % Compute source and destination regions for zoom mode:
            srcleft=min(xleft, xright);
            srctop=min(xtop, xbottom);
            srcright=max(xleft, xright);
            srcbottom=max(xtop, xbottom);
            
            if zoomset==2
                % Zoom set: Show selected region, expanded to full screen:
                texrect=CenterRect(Screen('Rect', texid), winrect);
                srcrect=[srcleft - texrect(1) srctop - texrect(2) srcright - texrect(1) srcbottom - texrect(2)];
                rscale = min(RectWidth(winrect) / RectWidth(srcrect), RectHeight(winrect) / RectHeight(srcrect));
                dstrect=CenterRect([0 0 RectWidth(srcrect)*rscale RectHeight(srcrect)*rscale], winrect);                
            else
                % No zoom: Show image 1-1:
                srcrect=[];
                dstrect=[];
                rscale=1;
            end
            
            % Draw selected srcrect region of our HDR image to selected
            % dstrect region of the screen:
            if hdrmode
                % Draw HDR texture:
                Screen('DrawTexture', win, texid, srcrect, dstrect, [], [], [], [sf sf sf]);
            else
                % Draw quantized 8bpc version:
                Screen('DrawTexture', win, ldrtexid, srcrect, dstrect, [], [], [], [sf sf sf]);
            end
            
            % Some status text:
            if hdrmode
                quantized = 'Floating point HDR';
            else
                quantized = 'Quantized to 8 bpc LDR';
            end
            
            DrawFormattedText(win, ['Image: ' imagename ' : Zoom=' num2str(rscale) ' : Intensity scaling: ' num2str(sf) ' : ' quantized ], 0, 0, [0 255 0]);

            % Show selected ROI during zoom selection, if any:
            if xright~=-1 && xbottom~=-1 && zoomset==1
                Screen('FrameRect', win, [255 0 0], [srcleft srctop srcright srcbottom]);
            end

            % Show updated HDR framebuffer at next vertical retrace:
            Screen('Flip', win);
                        
            % Clear backbuffer by overdrawing with a black full screen rect:
            Screen('FillRect', win, 0);

            % Keypress and what key?
            [pressed secs keycode] = KbCheck;
            if pressed
                % Key pressed: Dispatch...
                if keycode(KbName('space'))
                    % Load next image in folder:
                    needupdate = 0;

                    % Debounce key:
                    while KbCheck; end;
                end

                if keycode(KbName('q'))
                    % Switch between display of float32 HDR image and 8 bpc
                    % LDR image:
                    needupdate = 1;
                    hdrmode = 1 - hdrmode;
                    
                    % Debounce key:
                    while KbCheck; end;
                end
                
                if keycode(KbName('ESCAPE'))
                    % Exit viewer:
                    needupdate = 0;
                    abortit = 1;
                end

                if keycode(KbName('UpArrow'))
                    % Increase intensity scaling factor:
                    needupdate = 1;
                    sf = sf + sfstep;
                end

                if keycode(KbName('DownArrow'))
                    % Decrease intensity scaling factor:
                    needupdate = 1;
                    sf = sf - sfstep;
                end
            else
                % No key pressed: Check for mouse actions...
                [xm, ym, buttons]=GetMouse(win);

                if buttons(1)
                    % Left button pressed:
                    if zoomset==0
                        % Initial press: Reset and Setup zoom operation...
                        xleft=xm;
                        xtop=ym;
                        xright=-1;
                        xbottom=-1;
                        zoomset=1;
                    else
                        if zoomset==1
                            % Continued press: Update zoom area:
                            xright=xm;
                            xbottom=ym;
                        else
                            % Pressed while zoomed in: Reset...
                            while buttons(1)
                                [xm, ym, buttons]=GetMouse(win);
                            end
                            zoomset=0;                            
                        end
                    end
                else
                    % Left button released:
                    if zoomset==1
                        zoomset=2;
                    end
                end
            end
        end
        
        % Done with this image. Release its texture and next one...
        Screen('Close', texid);
        Screen('Close', ldrtexid);
        clear img;
        clear quantimg;
    end

    % Enable keypress output to Matlab window:
    ListenChar(0);
    
    % Release all textures, close all windows, shutdown BrightSide library:
    Screen('CloseAll');
    
    % Well done!
    fprintf('Bye bye!\n');
catch
    % Error handler: If something goes wrong between try and catch, we
    % close the window and abort.

    % Enable keypress output to Matlab window:
    ListenChar(0);

    % Release all textures, close all windows, shutdown BrightSide library:
    Screen('CloseAll');

    % Rethrow the error, so higher level routines can handle it:
    psychrethrow(psychlasterror);
end;
