function retval = SetAnaglyphStereoParameters(cmd, win, rgb, aux1, aux2)
% retval = SetAnaglyphStereoParameters(subcommand, windowPtr, [red green blue]);
%
% Change parameters of built-in anaglyph stereo display function at
% runtime. Allows to change all relevant settings of the algorithm at any
% time in your script. You need to open your onscreen window with the
% 'stereomode' parameter set to one of the anaglyph modes (i.e. one of
% 6=Red-Green anaglyph, 7=Green-Red, 8=Red-Blue, 9=Blue-Red, where
% Red-Green means "Left eye displayed in red channel, right eye displayed
% in green channel").
%
% Anaglyph stereo also works in stereo modes 2, 3, 4, 5 and 10 if you used
% the ...
% PsychImaging('AddTask', 'LeftView', 'DisplayColorCorrection', 'AnaglyphStereo')
% ... and ...
% PsychImaging('AddTask', 'RightView', 'DisplayColorCorrection', 'AnaglyphStereo')
% ... setup sequence to add an anaglyph shader to the end of each view
% channel. However, in those "non-pure" exotic anaglyph modes only a subset
% of the commands mentioned here work, specifically the ones for advanced
% anaglyph display, ie., subfunctions 'ColorAnaglyphMode' and friends. Inverted mode
% or simple mode or changing NTSC weights or simple gains does not work in
% this configuration.
%
% Additionally you also need to enable the Psychtoolbox image processing
% pipeline by using the PsychImaging() command for configuring and opening
% onscreen windows.
%
% Example:
% PsychImaging('PrepareConfiguration');
% stereomode = 6;
% windowPtr = PsychImaging('OpenWindow', screenid, bgcolor, [], [], [], stereomode);
%
% See 'help PsychImaging' and 'help PsychGLImageProcessing' for an overview
% of the imaging pipeline.
%
% Parameters and their meaning:
% -----------------------------
%
% 'windowPtr' is the handle to a anaglyph stereo onscreen window.
%
% 'subcommand' can be one of:
%
% 'LeftGains' - Provide (in optional vector [red green blue]) per-channel
% color gains for computation of the left-eye channel. This allows to
% adjust the output colors to optimally match the color emission / filter
% profile of your anaglyph glasses and monitor - Reduce ghosting / visual
% crosstalk. E.g.
%
% SetAnaglyphStereoParameters('LeftGains', win, [1.0 0.0 0.0]); sets the
% red-gain for the left eye image to 1.0 and disables output of the left
% eye image into the green or blue channels. This is the default setting
% for Red-Green or Red-Blue anaglyph, left eye view only goes into red
% channel.
%
% The command always returns the old gains as a 3 component vector.
%
% 'RightGains' - See left gains, but this time for right image:
% SetAnaglyphStereoParameters('RightGains', win, [0.0 0.7 0.0]);
% would set the green channel to 70% output intensity and disable the red-
% and blue channels. This is a typical setting for Red-Green anaglyph
% stereo: Right eye only written to green color channel, but with reduced
% output intensity to compensate for the stronger green-sensitivity of the
% human eye.
%
% 'ColorToLuminanceWeights' - Set the color to luminance conversion
% weights: PTB converts all color images into pure luminance (greyscale)
% images before distributing them into the color channels via the formula:
% Luminance = redvalue * redweight + bluevalue * blueweight + greenvalue *
% greenweight.
%
% The default red- green- and blue-weights are (0.3, 0.59, 0.11), so:
% L = 0.3*red + 0.59*green + 0.11*blue. This is according to NTSC standard
% (if i remember correctly!). You can change the weights via:
%
% SetAnaglyphStereoParameters('ColorToLuminanceWeights', win, [redweight greenweight blueweight]);
%
% The command always returns the old weights as a 3 component vector.
%
% 'BackgroundColorBias' - Set the "background color" for inverted anaglyph
% mode. This color value is added to each pixel. It defaults to zero - bias
% free display. For inverted mode, one usually sets it to maximum output
% intensity, so that the default background is a bright output.
%
% 'InvertedMode' - Switch to inverted display mode, i.e., set the
% 'BackgroundColorBias' and color gains properly. The image is displayed in
% an inverted fashion, like a photo-negative. This allows for less ghosting
% artifacts if your stimulus is sparse, e.g., dots and lines. This trick
% was proposed by Massimiliano di Luca.
%
% 'NormalMode' - Switch from inverted mode to normal anaglyph mode.
% Background color bias is set to zero, gains a positive.
%
% Color reproduction in color images with Anaglyphs:
%
% SetAnaglyphStereoParameters('ColorAnaglyphMode', win [,LeftMatrix] [,RightMatrix] [, Gamma]);
%
% This switches to a more flexible anaglyph implementation, where the
% output color is calculated as:
%
% [redout, greenout, blueout]' = LeftMatrix * [leftred, leftgreen, leftblue]' + RightMatrix * [rightred, rightgreen, rightblue]'
%
% where [redout, greenout, blueout] == Final anaglyph color pixel.
%       [leftred, leftgreen, leftblue] == Input left eye pixel color.
%       [rightred, rightgreen, rightblue] == Input right eye pixel color.
%
% and   LeftMatrix and RightMatrix are 3 by 3 color weight matrices.
%
% If 'Gamma' is provided, then the 'redout' color is gamma corrected with
% the given 'Gamma' value, ie. redout = redout ^ Gamma.
%
% Multiple sets of possible matrices and gamma weights are already defined
% and selectable with the following names:
%
% SetAnaglyphStereoParameters('OptimizedColorAnaglyphMode', win);
% SetAnaglyphStereoParameters('FullColorAnaglyphMode', win);
% SetAnaglyphStereoParameters('HalfColorAnaglyphMode', win);
% SetAnaglyphStereoParameters('GrayAnaglyphMode', win);
%
% The idea and matrices/gammas for these four selectable presets are taken
% from the following website:
%
% http://mitglied.lycos.de/stereo3d/anaglyphcomparison.htm
%
%
% Overlay support:
%
% If you want to have an additional monoscopic overlay window, e.g., for
% text instructions to the subject, fixation dots etc., you can get one via
% a call to:
%
% overlaywindow = SetAnaglyphStereoParameters('CreateMonoOverlay', win);
%
% This will create a fullscreen 'overlaywindow'. You can draw to that
% window like to any other window. It's content will be overlaid to the
% final anaglyph image, but without applying anaglyph conversion to it.
% Wherever this overlay has an alpha color value of zero, the stimulus will
% be shown. Wherever its alpha value is greater than zero, the overlays
% image will be shown.
%
% 'GetHandle' - Return GLSL handle to anaglyph shader. Allows to modify the
% shader itself, e.g., replace it by your own customized shader. Only for
% OpenGL experts!
%

% History:
% 5.1.2007   Wrote it (MK).
% 17.2.2007  Implemented setup code for Max di Lucas inversion trick (MK).
% 11.11.2007 Implemented setup code for extended/exotic anaglyph shading
%            and for color overlay support. (MK)
% 14.08.2012 Implement setup code for anaglyph stereo shaders inside per-view
%            channels, to handle separate display outputs, e.g., as on MPI
%            Kuka projection setup.

persistent inverted;

% Only works if GL shading language is supported.
if isempty(inverted)
    AssertGLSL;
    inverted = [];
end

try
    if nargin < 2
        error('Insufficient number of input arguments.');
    end
    
    % Assume combined shader in 'StereoCompositingBlit' chain:
    combinedAnaglyph = 1;

    % Query GLSL shader handle for anaglyph shader:
    [slot shaderid blittercfg voidptr glsl] = Screen('HookFunction', win, 'Query', 'StereoCompositingBlit', 'StereoCompositingShaderAnaglyph'); %#ok<*ASGLU>

    % Shader found?
    if slot == -1
        % Nope. Retry with per-view channels.
        combinedAnaglyph = 0;
        [slot shaderid blittercfg voidptr glsl(1)] = Screen('HookFunction', win, 'Query', 'StereoLeftCompositingBlit', 'ICM:AnaglyphStereo');
        if slot ~= -1
            [slot shaderid blittercfg voidptr glsl(2)] = Screen('HookFunction', win, 'Query', 'StereoRightCompositingBlit', 'ICM:AnaglyphStereo');
        end
    end
    
    % Still no shader?
    if slot == -1
        % Game over:
        error('Either the imaging pipeline is not enabled for given onscreen window, or it is not switched to Anaglyph stereo mode.');
    end
    
    if glsl == 0
        error('Anaglyph shader is not operational for unknown reason. Sorry...');
    end
    
    % Bind it:
    glUseProgram(glsl(1));
    
catch %#ok<*CTCH>
    % If anything failed. Unbind:
    glUseProgram(0);
    psychrethrow(psychlasterror);
end

    % Subcommand dispatch:
    if strcmpi(cmd, 'GetHandle')
        retval = glsl;
    end

    if strcmpi(cmd, 'InvertedMode')
        if ~combinedAnaglyph
            warning('PTB:SetAnaglyphStereoParametersUNSUPPORTED', 'SetAnaglyphStereoParameters(''InvertedMode'') is not supported for separate-viewchannel anaglyph mode. Ignored!');
            retval = 0;
            return;
        end
        
        % Switch to inverted display (Max di Lucas trick):
        if isempty(find(inverted == win)) %#ok<*EFIND>
            % Switch needed: Add windowhandle to our list of inverted
            % windows:
            inverted = [inverted win];
            
            % Retrieve background color bias, check if its non-zero:
            uniloc = glGetUniformLocation(glsl, 'ChannelBias');
            retval = glGetUniformfv(glsl, uniloc) * WhiteIndex(win);

            if max(abs(retval)) <= 0
                % Bias is zero, that's not suitable for inverted mode. Set
                % it to maximum output color, except for the components
                % where both gains are zero, i.e., unused channels:
                lg=SetAnaglyphStereoParameters('LeftGains', win);
                rg=SetAnaglyphStereoParameters('RightGains', win);
                glUseProgram(glsl);
                glUniform3fv(uniloc, 1, sign(lg + rg));
            end
            
            % Query and set color gains - This will switch gains properly:
            SetAnaglyphStereoParameters('LeftGains', win, SetAnaglyphStereoParameters('LeftGains', win));
            SetAnaglyphStereoParameters('RightGains', win, SetAnaglyphStereoParameters('RightGains', win));
            retval = 1;
        end
    end

    if strcmpi(cmd, 'StandardMode')
        if ~combinedAnaglyph
            warning('PTB:SetAnaglyphStereoParametersUNSUPPORTED', 'SetAnaglyphStereoParameters(''StandardMode'') is not supported for separate-viewchannel anaglyph mode. Ignored!');
            retval = 0;            
            return;
        end
        
        % Switch to standard, non-inverted display:
        winidx = find(inverted == win);
        if ~isempty(winidx)
            % Switch needed: Remove windowhandle from our list of inverted
            % windows, actually null it out, don't remove:
            inverted(winidx) = 0;
            
            % Retrieve background color bias, check if its non-zero:
            uniloc = glGetUniformLocation(glsl, 'ChannelBias');
            retval = glGetUniformfv(glsl, uniloc) * WhiteIndex(win);
            
            if max(abs(retval)) > 0
                % Bias is non-zero, that's not suitable for inverted mode. Set
                % it to black output color:
                glUniform3fv(uniloc, 1, [0 0 0]);
            end
            
            % Query and set color gains - This will switch gains properly:
            SetAnaglyphStereoParameters('LeftGains', win, SetAnaglyphStereoParameters('LeftGains', win));
            SetAnaglyphStereoParameters('RightGains', win, SetAnaglyphStereoParameters('RightGains', win));
            retval = 0;
        end
    end

    if strcmpi(cmd, 'BackgroundColorBias')
        if ~combinedAnaglyph
            warning('PTB:SetAnaglyphStereoParametersUNSUPPORTED', 'SetAnaglyphStereoParameters(''BackgroundColorBias'') is not supported for separate-viewchannel anaglyph mode. Ignored!');
            retval = [0, 0, 0];
            return;
        end
        
        uniloc = glGetUniformLocation(glsl, 'ChannelBias');
        retval = glGetUniformfv(glsl, uniloc) * WhiteIndex(win);
        if nargin>=3
            if size(rgb)~=3
                error('Provided call parameter must be a 3 component vector with (R,G,B) bias color values.');
            end
            
            % Normalize from color range of window to GL's 0-1 range:
            rgb = rgb / WhiteIndex(win);
            
            glUniform3fv(uniloc, 1, rgb);
        end
    end

    if strcmpi(cmd, 'ColorToLuminanceWeights')
        if ~combinedAnaglyph
            warning('PTB:SetAnaglyphStereoParametersUNSUPPORTED', 'SetAnaglyphStereoParameters(''ColorToLuminanceWeights'') is not supported for separate-viewchannel anaglyph mode. Ignored!');
            retval = [0, 0, 0];
            return;
        end
        
        uniloc = glGetUniformLocation(glsl, 'ColorToGrayWeights');
        retval = glGetUniformfv(glsl, uniloc);
        if nargin>=3
            if size(rgb)~=3
                error('Provided call parameter must be a 3 component vector with color weights or gains.');
            end
            glUniform3fv(uniloc, 1, rgb);
        end
    end

    if strcmpi(cmd, 'LeftGains')
        if ~combinedAnaglyph
            warning('PTB:SetAnaglyphStereoParametersUNSUPPORTED', 'SetAnaglyphStereoParameters(''LeftGains'') is not supported for separate-viewchannel anaglyph mode. Ignored!');
            retval = [0, 0, 0];
            return;
        end
        
        uniloc = glGetUniformLocation(glsl, 'Gains1');
        retval = abs(glGetUniformfv(glsl, uniloc));

        if find(inverted ==  win)
            inverter = -1;
        else
            inverter = 1;
        end

        if nargin>=3
            if size(rgb)~=3
                error('Provided call parameter must be a 3 component vector with color weights or gains.');
            end
            glUniform3fv(uniloc, 1, inverter * rgb);
        end
    end
    
    if strcmpi(cmd, 'RightGains')
        if ~combinedAnaglyph
            warning('PTB:SetAnaglyphStereoParametersUNSUPPORTED', 'SetAnaglyphStereoParameters(''RightGains'') is not supported for separate-viewchannel anaglyph mode. Ignored!');
            retval = [0, 0, 0];
            return;
        end
        
        uniloc = glGetUniformLocation(glsl, 'Gains2');
        retval = abs(glGetUniformfv(glsl, uniloc));

        if find(inverted == win)
            inverter = -1;
        else
            inverter = 1;
        end

        if nargin>=3
            if size(rgb)~=3
                error('Provided call parameter must be a 3 component vector with color weights or gains.');
            end
            glUniform3fv(uniloc, 1, inverter * rgb);
        end
    end

    if strcmpi(cmd, 'FullColorAnaglyphMode')
        % Implement colored anaglyph rendering for full color anaglyphs:
        
        % Store left-gain matrix in 'rgb':
        rgb =  [[1 0 0]; [0 0 0]; [0 0 0]];
        
        % Store right-gain matrix in 'aux1':
        aux1 = [[0 0 0]; [0 1 0]; [0 0 1]];
        
        % Rewrite cmd, so the generic setup code for color anaglyphs gets
        % called to do the actual setup work:
        cmd = 'ColorAnaglyphMode';
    end

    if strcmpi(cmd, 'GrayAnaglyphMode')
        % Implement colored anaglyph rendering for gray anaglyphs:
        
        % Store left-gain matrix in 'rgb':
        rgb =  [[0.299 0.587 0.114]; [0 0 0]; [0 0 0]];
        
        % Store right-gain matrix in 'aux1':
        aux1 = [[0 0 0]; [0.299 0.587 0.114]; [0.299 0.587 0.114]];
        
        % Rewrite cmd, so the generic setup code for color anaglyphs gets
        % called to do the actual setup work:
        cmd = 'ColorAnaglyphMode';
    end

    if strcmpi(cmd, 'HalfColorAnaglyphMode')
        % Implement colored anaglyph rendering for half color anaglyphs:
        
        % Store left-gain matrix in 'rgb':
        rgb =  [[0.299 0.587 0.114]; [0 0 0]; [0 0 0]];
        
        % Store right-gain matrix in 'aux1':
        aux1 = [[0 0 0]; [0 1 0]; [0 0 1]];
        
        % Rewrite cmd, so the generic setup code for color anaglyphs gets
        % called to do the actual setup work:
        cmd = 'ColorAnaglyphMode';
    end

    if strcmpi(cmd, 'OptimizedColorAnaglyphMode')
        % Implement colored anaglyph rendering for optimized color anaglyphs:
        
        % Store left-gain matrix in 'rgb':
        rgb =  [[0 0.7 0.3]; [0 0 0]; [0 0 0]];
        
        % Store right-gain matrix in 'aux1':
        aux1 = [[0 0 0]; [0 1 0]; [0 0 1]];
        
        % Setup gamma correction for red channel with a gamma of 1.5:
        aux2 = 1.5;
        
        % Rewrite cmd, so the generic setup code for color anaglyphs gets
        % called to do the actual setup work:
        cmd = 'ColorAnaglyphMode';
    end

    if strcmpi(cmd, 'ColorAnaglyphMode')
        % Switch from standard luminance monochrome anaglyph shader (as
        % created by Screen('OpenWindow', ...) by default) to a more
        % complex and flexible shader. Either that, or reparameterize that
        % shader...
        
        % Color anaglyph shader already bound to imaging pipeline? Only for
        % regular combined case, as we know the proper shader is bound
        % otherwise, because it was already created by
        % PsychColorCorrection() during pipeline setup in PsychImaging():
        if combinedAnaglyph && ~strcmpi(shaderid, 'StereoCompositingShaderAnaglyphExtended')
            % Nope. Destroy old default shader and replace it by our
            % shader:
            Screen('HookFunction', win, 'Remove', 'StereoCompositingBlit', slot);
            
            % Load shader from file:
            colorshader = LoadGLSLProgramFromFiles('ColoredAnaglyphsShader', 1);
            glsl = colorshader;

            % Setup stereo image to texture unit mappings:
            glUseProgram(colorshader);
            glUniform1i(glGetUniformLocation(colorshader, 'Image1'), 0);
            glUniform1i(glGetUniformLocation(colorshader, 'Image2'), 1);

            % Setup default gamma correction for red channel to zero, ie,
            % no gamma correction at all!
            glUniform1f(glGetUniformLocation(colorshader, 'RedGamma'), 0);
            
            % Insert it at former position of the old shader:
            posstring = sprintf('InsertAt%iShader', slot);
            Screen('Hookfunction', win, posstring, 'StereoCompositingBlit', 'StereoCompositingShaderAnaglyphExtended', colorshader);
            
            % Ok, shader ready. Now for the setup of conversion
            % parameters...
        end
        
        if combinedAnaglyph
            % Regular case, one shader to rule them all:
            unilocleft = glGetUniformLocation(glsl, 'GainsLeft');
            unilocright = glGetUniformLocation(glsl, 'GainsRight');
            gammaloc = glGetUniformLocation(glsl, 'RedGamma');
            
            if exist('rgb', 'var') && ~isempty(rgb)
                if size(rgb,1)~=3 || size(rgb,2)~=3
                    error('Provided left-gain matrix parameter must be a 3x3 component matrix.');
                end
                glUniformMatrix3fv( unilocleft, 1, 0, rgb);
            end
            
            if exist('aux1', 'var') && ~isempty(aux1)
                if size(aux1,1)~=3 || size(aux1,2)~=3
                    error('Provided right-gain matrix parameter must be a 3x3 component matrix.');
                end
                glUniformMatrix3fv( unilocright, 1, 0, aux1);
            end
            
            if exist('aux2', 'var') && ~isempty(aux2)
                if size(aux2,1)~=1 || size(aux2,2)~=1
                    error('Provided gamma parameter must be a scalar.');
                end
                glUniform1f( gammaloc, aux2);
            end
        else
            % Separate shaders for left- and right view channel:

            % Handle left channel shader first:
            unilocleft = glGetUniformLocation(glsl(1), 'GainsLeft');
            gammaloc = glGetUniformLocation(glsl(1), 'RedGamma');
            
            glUseProgram(glsl(1));
            if exist('rgb', 'var') && ~isempty(rgb)
                if size(rgb,1)~=3 || size(rgb,2)~=3
                    error('Provided left-gain matrix parameter must be a 3x3 component matrix.');
                end
                glUniformMatrix3fv( unilocleft, 1, 0, rgb);
            end
            
            if exist('aux2', 'var') && ~isempty(aux2)
                if size(aux2,1)~=1 || size(aux2,2)~=1
                    error('Provided gamma parameter must be a scalar.');
                end
                glUniform1f( gammaloc, aux2);
            end

            % Handle right channel shader similarly:
            glUseProgram(glsl(2));

            % Even the "right gain matrix" is named "GainsLeft" in this case:
            unilocright = glGetUniformLocation(glsl(2), 'GainsLeft');
            gammaloc = glGetUniformLocation(glsl(2), 'RedGamma');
            
            if exist('aux1', 'var') && ~isempty(aux1)
                if size(aux1,1)~=3 || size(aux1,2)~=3
                    error('Provided right-gain matrix parameter must be a 3x3 component matrix.');
                end
                glUniformMatrix3fv( unilocright, 1, 0, aux1);
            end
            
            if exist('aux2', 'var') && ~isempty(aux2)
                if size(aux2,1)~=1 || size(aux2,2)~=1
                    error('Provided gamma parameter must be a scalar.');
                end
                glUniform1f( gammaloc, aux2);
            end            
        end
    end

    if strcmpi(cmd, 'CreateMonoOverlay')
        winfo = Screen('GetWindowInfo', win);
        if winfo.StereoMode < 2 || winfo.StereoMode > 9
            error('SetAnaglyphStereoParameters(''CreateMonoOverlay'') is not supported in current stereomode %i!', winfo.StereoMode);
        end

        % Create a monoscopic overlay window which is not affected by
        % anaglyph conversion, but just overlaid to the final rendering.
        glUseProgram(0);
        
        % Create Offscreen window for the overlay. It has the same size as
        % the onscreen window, the same pixeldepth, but a completely black
        % background with alpha value zero -- fully transparent by default.
        overlaywin = Screen('OpenOffscreenWindow', win, [0 0 0 0]);
        
        % Retrieve low-level OpenGl texture handle to the window:
        overlaytex = Screen('GetOpenGLTexture', win, overlaywin);

        % Append blitter command for a one-to-one blit of the overlay
        % texture to the target buffer. We need to disable the 2nd texture
        % unit and setup proper alpha testing, so a) it works at all and b)
        % the overly only occludes the main image where the overly has a
        % non-zero alpha:
        Screen('Hookfunction', win, 'AppendMFunction', 'StereoCompositingBlit', 'Setup1 Alphatest for Overlay', 'glAlphaFunc(516, 0.0);');
        Screen('Hookfunction', win, 'AppendMFunction', 'StereoCompositingBlit', 'Setup2 Alphatest for Overlay', 'glEnable(3008);');
        Screen('Hookfunction', win, 'AppendMFunction', 'StereoCompositingBlit', 'Setup3 Texunit1 off for Overlay', 'glActiveTexture(33985);');
        Screen('Hookfunction', win, 'AppendMFunction', 'StereoCompositingBlit', 'Setup4 Texunit1 off for Overlay', 'glDisable(34037);');
        Screen('Hookfunction', win, 'AppendMFunction', 'StereoCompositingBlit', 'Setup5 Texunit1 off for Overlay', 'glActiveTexture(33984);');
        Screen('Hookfunction', win, 'AppendBuiltin',   'StereoCompositingBlit', 'Builtin:IdentityBlit', sprintf('TEXTURERECT2D(0)=%i', overlaytex));
        Screen('Hookfunction', win, 'AppendMFunction', 'StereoCompositingBlit', 'Teardown Alphatest for Overlay', 'glDisable(3008);');
        
        % Return handle to overlay window:
        retval = overlaywin;
    end
    
try
    % Unbind shader again...
    glUseProgram(0);
catch
    psychrethrow(psychlasterror);
end

% Done.
return;
