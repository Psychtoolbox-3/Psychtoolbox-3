function DeinterlacerTest(sf, imgname, method)
% DeinterlacerTest([sf][,imgname][,method])
%
% Tests (and thereby demonstrates) the Psychtoolbox builtin GLSL realtime
% deinterlaceshader, an OpenGL GLSL Fragmentshader for on-the-fly
% deinterlacing of video images which are given as PTB textures and drawn
% via Screen('MakeTexture').
%
% Parameters: 'imgname' is optional. If provided, the routine loads an
% imagefile of that name as test input, otherwise it creates a synthetic
% test pattern in Matlab.
%
% sf - Optional scalefactor, by which the image is scaled in size.
%
% method - Select deinterlace method:
%
% 0 = None. 1 = Replacement of lines, 2 = Averaging of lines.
%
% See the shader sources under Psychtoolbox/PsychOpenGL/PsychGLSLShaders
% for proper definition of the operations.
%
% What you should see (assuming you did not provide an input image file):
%
% First an image with alternating red- and green color gradients. This is
% the interlaced test pattern where even rows contain only red, odd rows
% contain only green gradients.
%
% Press a key.
%
% Then, by pressing and releasing a key, the image should alternate between
% an all red image (only red gradients of twice the thickness), and an all
% green image: This is the alternating deinterlaced even, odd, even, odd
% ... half-fields of the interlaced image.
%
% End the demo by pressing ESCape key twice.
%
% If you see something else, or just a black screen, or it aborts with
% error, then your hardware is probably not capable of running the
% deinterlace-shader and you're out of luck.
%

AssertOpenGL;
filter=0;

KbName('UnifyKeyNames');
escape=KbName('ESCAPE');

try
    screenid = max(Screen('Screens'));
    [win winrect] = Screen('OpenWindow', screenid, 0);
    
    if nargin<3
        method=1;
    end
    
    if nargin<1
        sf = 1;
    end
    
    if nargin<2
        imgname=[];
    end
    
    % Image filename provided as input?
    if isempty(imgname)
        % Nope. Create synthetic test image:
        img = uint8(zeros(256,256,3));
        for y=0:size(img,1)-1
            for x=1:256
                % Even row: Red gradient
                img(y+1, x, mod(y,2) + 1)=(x-1);
            end
            
            if mod(floor(y/2),2)==0
                img(y+1, :, mod(y,2) + 1)=255 - img(y+1, :, mod(y,2) + 1);
            end
        end
    else
        % Read image from filesystem:
        img = imread(imgname);
    end
    
    for i=1:3
        teximg(:,:,i)=transpose(img(:,:,i));
    end
    
    % Build corresponding texture:
    tex = Screen('MakeTexture', win, teximg);
    
    dstrect = Screen('Rect', tex);
    dstrect = ScaleRect(dstrect, sf, sf);
    dstrect = CenterRect(dstrect, winrect);
    
    % Show original (interlaced) image:
    Screen('DrawTexture', win, tex, [], dstrect, 90, filter);
    Screen('Flip', win);
    
    % Wait for keypress:
    KbStrokeWait;
    
    % Load deinterlace-shader:
    deinterlacer = [];
    if method == 1
        % Simple line replication:
        deinterlacer = LoadGLSLProgramFromFiles('EXPDeinterlaceShaderLineDouble',1);
    else
        if method ==2
            % Replace line by average of neighbour lines:
            deinterlacer = LoadGLSLProgramFromFiles('EXPDeinterlaceShaderAvg',1);
        end
    end
    
    if method > 0
        % Bind and initialize it:
        glUseProgram(deinterlacer);
        % Input image will be bound to texture unit zero:
        glUniform1i(glGetUniformLocation(deinterlacer, 'Image1'), 0);
        % Get handle for the field selection parameter:
        useoddfield=glGetUniformLocation(deinterlacer, 'UseOddField');
    end
    
    % Run until user presses ESCape or for 1000 iterations.
    count = 0;
    while count < 1000
        % Show the even half-field: Tell deinterlacer we want the even field:
        if method>0
            glUseProgram(deinterlacer);
            glUniform1f(useoddfield, 0);
        end
        
        % Draw it:
        Screen('DrawTexture', win, tex, [], dstrect, 90, filter, [], [], deinterlacer);
        Screen('Flip', win, 0, 0, 0);
        
        % Wait for keypress:
        KbStrokeWait;
        
        % Show the odd half-field: Tell deinterlacer we want the odd field:
        if method>0
            glUseProgram(deinterlacer);
            glUniform1f(useoddfield, 1);
        end
        
        % Draw it:
        Screen('DrawTexture', win, tex, [], dstrect, 90, filter, [], [], deinterlacer);
        Screen('Flip', win, 0, 0, 0);
        
        % Wait for keypress:
        [dummy keycode] = KbStrokeWait;
        if keycode(escape)
            break;
        end
        
        count = count + 1;
    end
    
    % Shut down:
    Screen('CloseAll');
catch
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
