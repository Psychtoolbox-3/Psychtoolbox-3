function [warpstruct, filterMode] = CreateDisplayWarp(window, calibfilename, showCalibOutput, varargin)
% [warpstruct, filterMode] = CreateDisplayWarp(window, calibfilename [, showCalibOutput=0]);
%
% Helper routine for Geometric display undistortion mapping, not to be
% called inside normal PTB scripts!
%
% This function reads a display calibration file 'calibfilename' and builds
% a "geometric warp function" based on the calibration information in
% 'calibfilename' for the onscreen window with handle 'window'. It returns
% a struct 'warpstruct' that defines the created warp function. You could
% pass this 'warpstruct' as a parameter to the Psychtoolbox command...
%
% PsychImaging('AddTask', viewchannel, 'GeometryCorrection', warpstruct);
%
% However, you normally do not call this routine directly from your script. Its
% called internally by the PsychImaging() command...
% 
% PsychImaging('AddTask', viewchannel, 'GeometryCorrection', calibfilename);
%
% ...in order to setup PTB's imaging pipeline for realtime geometry
% correction, based on the calibration info in the file 'calibfilename'.
%
% Example: You created a calibration file 'mycalib.mat' to undistort the
% left view display of a stereo setup. Then you could apply this
% undistortion function via the following setup code:
%
% PsychImaging('PrepareConfiguration');
% PsychImaging('AddTask', 'LeftView', 'GeometryCorrection', 'mycalib.mat');
% window = PsychImaging('OpenWindow', screenid);
%
% This would open an onscreen window just as window=Screen('OpenWindow', screenid);
% would do. It would configure the window for automatic undistortion based
% on the data in 'mycalib.mat'.
%
% Psychtoolbox provides multiple different interactive setup routines that
% allow you to create a calibration file for your setup. Currently the
% following routines are provided:
%
% DisplayUndistortionBVL.m    -- Undistortion based on 3rd order polynomial
% surface. This is the recommended calibration procedure for most cases -
% Proven in real-world use on many different display types.
%
% DisplayUndistortionBezier.m -- Undistortion based on a NURBS surface (Non
% uniform rational bezier spline surface). A simple procedure.
%
% DisplayUndistortionHalfCylinder.m -- Undistortion for projection of
% images to a half-cylindrical projection surface.
%
% DisplayUndistortionSphere.m -- Undistortion for projection of
% images to a spherical or half-spherical projection surface.
%
% DisplayUndistortionCSV.m -- Import undistortion information from
% a .csv file with a warp mesh description suitable for use with NVidia's
% Warp API. This creates a compatible display warping to use of that NVidia
% technology.
%

% History:
% 19.7.2007 Written (MK).
% 17.2.2008 Added undistortion method donated by the Banks Vision Lab (MK).
% 10.3.2008 Fixed image inversion bug in BVL calibration (MK).
%  2.5.2008 Add support for bilinear texture filter shader to handle float
%           framebuffers on hw that doesn't filter float textures (MK).
% 13.4.2009 Improved support for bilinear texture filter shaders. (MK).
%           Optional 'Query' command to query last warpstruct.    (MK).
%           Support for half-cylinder projection. (MK).
% 25.8.2011 Adapt code for sphere projection undistortion to new convention
%           of Ingmar Schneider's shader code. (MK).
% 27.7.2012 Add support for DisplayUndistortionCSV() aka "NVidia Warp-API" format (MK).
%

% Global GL handle for access to OpenGL constants needed in setup:
global GL;

% Cache last generated warpstruct, so code can easily query it:
persistent oldwarpstruct;

if isempty(GL)
    sca;
    error('PTB internal error: GL struct not initialized?!?');
end

% Special case of simple query of last created 'warpstruct'?
if nargin == 1
    if ~ischar(window)
        error('Single provided argument is not a command string!');
    end
    
    if ~strcmpi(window, 'Query')
        error('Single provided argument is not the command string ''Query''!');
    end
    
    % "Query" command recognized. Return last created warpstruct:
    warpstruct = oldwarpstruct;
    return;
end

if nargin < 2
    sca;
    error('PTB internal error: Must provide all parameters!');
end

if nargin < 3 || isempty(showCalibOutput)
    showCalibOutput = 0;
end

% Is calibfilename a struct with calibration settings, or a filename of a
% calibration file?
if isstruct(calibfilename)
    % A struct: Assign it directly.
    calib = calibfilename;
else
    % Supposedly the filename of a calibration file:
    if ~ischar(calibfilename)
        error('In setup of geometry undistortion: Parameter "calibfilename" is not a filename string!');
    end
    
    % Load calibration file:
    if ~exist(calibfilename, 'file')
        sca;
        error('In setup of geometry undistortion: No such calibration file %s!', calibfilename);
    end

    calib = load(calibfilename);
end

% Preinit warpstruct:
warpstruct.glsl = [];
warpstruct.gld = [];

% Assume no need for texture filter shader:
needFilterShader = 0;
filterMode = ':Bilinear';

% Do we need a GLSL texture filter shader? We'd need one if the given
% gfx-hardware is not capable of filtering the input image buffer:
winfo = Screen('GetWindowInfo', window);
effectivebpc = 8;
if winfo.BitsPerColorComponent >= 16
    % Window is a floating point window with at least 16bpc.
    effectivebpc = 16;
    
    if winfo.BitsPerColorComponent >= 32
        % All buffers are 32 bpc for certain:
        effectivebpc = 32;
    end
    
    if (winfo.BitsPerColorComponent == 16)
        % First buffer is 16 bpc, following ones could be 32 bpc:
        if bitand(winfo.ImagingMode, kPsychUse32BPCFloatAsap)
            % All following buffers are 32bpc float. In the tradition of
            % "better safe than sorry", we assume that the warp op will use
            % one of the 32 bpc float buffers as input.
            effectivebpc = 32;            
        end
    end    
end

% Highres input buffer?
if effectivebpc > 8
    % Yes. Our input is a float texture. Check if the hardware can filter
    % textures of effectivebpc bpc in hardware:
    if effectivebpc > winfo.GLSupportsFilteringUpToBpc
        % Hardware not capable of handling such deep textures. We need to
        % create and attach our own bilinear texture filter shader:
        needFilterShader = 1;
        filterMode = '';
    end
end

% Actual setup code for display warp struct.
% ==========================================

% Type of setup depends on type of calibration:
switch(calib.warptype)
    case {'HalfCylinderProjection', 'SphereProjection'}
        % Build combo of displaylist and GLSL shader for projection of flat
        % screen image onto a half-cylinder or sphere:

        % Query effective onscreen window size:
        [winWidth, winHeight] = Screen('WindowSize', window);
                
        % Build the unwarp mesh display list within the OpenGL context of
        % Psychtoolbox:
        Screen('BeginOpenGL', window, 1);
                
        % Build a display list that corresponds to the current calibration:
        gld = glGenLists(1);
        glNewList(gld, GL.COMPILE);
        
        % "Draw" the warp-mesh once, so it gets recorded in the display list:
        if isempty(calib.rotationAngle)
            calib.rotationAngle = 0;
        end
        
        if isempty(calib.inSize)
            calib.inSize = [winWidth, winHeight];
        end
        
        if isempty(calib.inOffset)
            calib.inOffset = [0, 0];
        end
        
        if isempty(calib.outOffset)
            calib.outOffset = [0, 0];
        end
        
        if isempty(calib.outSize)
            calib.outSize = [winWidth, winHeight];
        end

        if isempty(calib.Wflat)
            calib.Wflat = 44;
        end

        if isempty(calib.R)
            calib.R = 32;
        end

        % No color gain correction:
        glColor4f(1,1,1,1);
        
        glTranslatef(calib.outOffset(1), calib.outOffset(2), 0);
        
        % Apply some rotation correction for misaligned displays:
        glTranslatef(calib.outSize(1)/2, calib.outSize(2)/2, 0);
        glRotatef(calib.rotationAngle, 0.0, 0.0, 1.0);
        glTranslatef(-calib.outSize(1)/2, -calib.outSize(2)/2, 0);
        
        % Draw a single default quad:
        glBegin(GL.QUADS)
        glTexCoord2f(0,calib.outSize(2));
        glVertex2f(0,0);

        glTexCoord2f(calib.outSize(1),calib.outSize(2));
        glVertex2f(calib.outSize(1),0);

        glTexCoord2f(calib.outSize(1),0);
        glVertex2f(calib.outSize(1),calib.outSize(2));

        glTexCoord2f(0,0);
        glVertex2f(0,calib.outSize(2));
        glEnd;

        
        % List ready - and already updated in the imaging pipeline:
        glEndList;

        Screen('EndOpenGL', window);

        % Assign display list to output warpstruct:
        warpstruct.gld = gld;

        if strcmpi(calib.warptype, 'SphereProjection')
            % Use spherical projection shader:
            warpstruct.glsl = LoadGLSLProgramFromFiles('SphereProjectionShader');
        else
            % Use cylindrical projection shader:
            warpstruct.glsl = LoadGLSLProgramFromFiles('CylinderProjectionShader');
        end
        glUseProgram(warpstruct.glsl);
        
        glUniform1i(glGetUniformLocation(warpstruct.glsl, 'doFilter'), needFilterShader);
        glUniform1i(glGetUniformLocation(warpstruct.glsl, 'Image'), 0);
        glUniform2f(glGetUniformLocation(warpstruct.glsl, 'inSize'), calib.inSize(1), calib.inSize(2));
        glUniform2f(glGetUniformLocation(warpstruct.glsl, 'inOffset'), calib.inOffset(1), calib.inOffset(2));
        glUniform2f(glGetUniformLocation(warpstruct.glsl, 'outSize'), calib.outSize(1), calib.outSize(2));

        if strcmpi(calib.warptype, 'SphereProjection')
            % Additional parameters for sphere projection:
            glUniform1f(glGetUniformLocation(warpstruct.glsl, 'Wflat'), calib.Wflat);
            glUniform1f(glGetUniformLocation(warpstruct.glsl, 'R'), calib.R);
        end        
        
        glUseProgram(0);
        
        % Ready.
        
    case {'BezierDisplayList'}
        % Build warp display list for Bezier surface based
        % calibration/remapping:
        
        Screen('BeginOpenGL', window, 1);
        
        gld = glGenLists(1);
        
        % Build a display list that corresponds to the current calibration:
        glNewList(gld, GL.COMPILE);

        glColor4f(1,1,1,1);

        subdivision = calib.subdivision;
        
        % Setup a 2D parametric grid with 'subdivision' subdivisions:
        glMapGrid2d(subdivision, 0, 1, subdivision, 0, 1);

        % Enable Bezier evaluators:
        glEnable(GL.MAP2_VERTEX_3);
        glEnable(GL.MAP2_TEXTURE_COORD_2);

        % Setup initial mapping table for texture coordinates (source image control
        % points):
        frompts = calib.frompts;
        % Establish mapping for texture coordinates:
        glMap2d(GL.MAP2_TEXTURE_COORD_2, 0, 1, 2, size(frompts,2), 0, 1, 2*size(frompts,2), size(frompts,3), frompts);

        % Establish mapping for vertex coordinates:
        topts = calib.topts;
        % Setup mapping based on current control point matrix for destination
        % points:
        glMap2d(GL.MAP2_VERTEX_3, 0, 1, 3, size(topts,2), 0, 1, 3*size(topts,2), size(topts,3), topts);
        
        % Compute the mesh based on current mappings:
        glEvalMesh2(GL.FILL, 0, subdivision, 0, subdivision);
        
        % Disable mesh evaluators:
        glDisable(GL.MAP2_VERTEX_3);
        glDisable(GL.MAP2_TEXTURE_COORD_2);

        % List ready - and already updated in the imaging pipeline:
        glEndList;

        % Assign display list to output warpstruct:
        warpstruct.gld = gld;

        % Ready.
        Screen('EndOpenGL', window);
        
    case {'BVLDisplayList'}
        % Build warp display list for calibration/remapping method
        % donated by the Banks Vision Lab:
        
        % Query effective onscreen window size:
        [winWidth, winHeight] = Screen('WindowSize', window);
        
        % Additional optional parameters provided?
        xLoomSize = [];
        yLoomSize = [];

        if nargin >= 5
            % At least two additional parameters. Really?
            if ~isempty(varargin{1})
                xLoomSize = varargin{1};
            end
            
            if ~isempty(varargin{2})
                yLoomSize = varargin{2};
            end
        end
        
        % Compute vertex- and texcoord-arrays that define the mesh
        % of quadrilaterals which should be rendered (with the stimulus
        % texture applied) to create the undistortion warp:
        [xyzcalibpos, xytexcoords] = BVLComputeWarpMesh(winWidth, winHeight, calib.scal, showCalibOutput, xLoomSize, yLoomSize);

        % Build the unwarp mesh display list within the OpenGL context of
        % Psychtoolbox:
        Screen('BeginOpenGL', window, 1);
                
        % Build a display list that corresponds to the current calibration:
        gld = glGenLists(1);
        glNewList(gld, GL.COMPILE);
        
        % "Draw" the warp-mesh once, so it gets recorded in the display list:
        glColor4f(1,1,1,1);
        glEnableClientState(GL.VERTEX_ARRAY);
        glVertexPointer(2, GL.DOUBLE, 0, xyzcalibpos);
        glEnableClientState(GL.TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL.DOUBLE, 0, xytexcoords);

        glDrawArrays(GL.QUADS, 0, length(xyzcalibpos)/2);

        glDisableClientState(GL.TEXTURE_COORD_ARRAY);
        glDisableClientState(GL.VERTEX_ARRAY);

        % List ready - and already updated in the imaging pipeline:
        glEndList;

        Screen('EndOpenGL', window);

        % Assign display list to output warpstruct:
        warpstruct.gld = gld;

        % Ready.
        
    case {'CSVDisplayList'}
        % Build warp display list for a calibration/remapping method that
        % is compatible with the method used, e.g., by NVidia's Warp API:
        
        % Query effective onscreen window size:
        [winWidth, winHeight] = Screen('WindowSize', window);
        
        % Compute vertex- and texcoord-arrays that define the mesh
        % of quadrilaterals which should be rendered (with the stimulus
        % texture applied) to create the undistortion warp:
        [xyzcalibpos, xytexcoords] = CSVComputeWarpMesh(winWidth, winHeight, calib.scal, showCalibOutput);

        % Build the unwarp mesh display list within the OpenGL context of
        % Psychtoolbox:
        Screen('BeginOpenGL', window, 1);
        
        % Build a display list that corresponds to the current calibration:
        gld = glGenLists(1);
        glNewList(gld, GL.COMPILE);
        
        % "Draw" the warp-mesh once, so it gets recorded in the display list:
        glColor4f(1,1,1,1);
        glEnableClientState(GL.VERTEX_ARRAY);
        glVertexPointer(2, GL.DOUBLE, 0, xyzcalibpos);
        glEnableClientState(GL.TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL.DOUBLE, 0, xytexcoords);

        glDrawArrays(GL.QUADS, 0, length(xyzcalibpos)/2);

        glDisableClientState(GL.TEXTURE_COORD_ARRAY);
        glDisableClientState(GL.VERTEX_ARRAY);

        % List ready - and already updated in the imaging pipeline:
        glEndList;

        Screen('EndOpenGL', window);

        % Assign display list to output warpstruct:
        warpstruct.gld = gld;

        % Ready.
        
    otherwise
        sca;
        error('Unknown calibration method id: %s!', calib.warptype);
end

% Need a filtershader and don't have one assigned yet?
if (needFilterShader > 0) && isempty(warpstruct.glsl)
    % Yes. Load, create and assign our default bilinear texture filter
    % shader:
    warpstruct.glsl = LoadGLSLProgramFromFiles('BilinearTextureFilterShader');
    glUseProgram(warpstruct.glsl);
    glUniform1i(glGetUniformLocation(warpstruct.glsl, 'Image'), 0);
    glUseProgram(0);
end

% Cache created warptstruct for later queries:
oldwarpstruct = warpstruct;

% Done. Return the warpstruct:
return;

% --- Helper routines for setup of the calibration method 'BVLDisplayList' ---

function [xyzcalibpos, xytexcoords] = BVLComputeWarpMesh(windowWidth, windowHeight, scal, showCalibOutput, xLoomSize, yLoomSize)
% [xyzcalibpos, xytexcoords] = BVLComputeWarpMesh(windowWidth, windowHeight, scal, showCalibOutput, xLoomSize, yLoomSize)
%
% Internal helper routine: Called by CreateDisplayWarp.m, which is in turn
% called by PsychImaging.m. Implements the geometric display calibration
% and undistortion procedure developed by the Banks Vision Lab at UC
% Berkeley.
%
% Use the calibration information stored in 'scal', together with the
% current 'windowWidth' and 'windowHeight' of the onscreen window to
% calibrate and compute vectors of vertex coordinates and texture
% coordinates for a mesh that performs the proper "display undistortion".
%
% History:
% 02/17/08  Derived from BFloadtimecalib.m with minimal modifications. (MK)
%

% Check resolution against the calibration file resolutions:
if (windowWidth ~= scal.rect(3)) || (windowHeight ~= scal.rect(4))
    fprintf('\n\nCALIBRATION WARNING!\n');
    fprintf('Onscreen window resolution (%d, %d) does not match ', ...
        windowWidth, windowHeight);
    fprintf('the resolution used in the calibration file (%d, %d)!\n', ...
        scal.rect(3), scal.rect(4));
    fprintf('Using the window resolution to draw the stimuli.\n');
    fprintf('\n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration
% Generate the calibration vertices here so we only do this once.
% We calculate the uncalibrated vertex coordinates, in screenspace
% (pixel coordinates), then pass them into our fitting routine to
% interpolate the calibrated coordinates.

% PreCompute vertex values
if isempty(xLoomSize)
    xLoomSize   = 73; %length(XVALUES)   %Can reduce loom resolution to speed up calib
end

if isempty(yLoomSize)
    yLoomSize   = 53; %length(YVALUES)
end

xStep       = windowWidth / (xLoomSize-1);
yStep       = windowHeight / (yLoomSize-1);

numVerts    = xLoomSize * yLoomSize;
% vertexCoords = Nx2 array, N rows of [x y] pairs.
vertexCoords    = zeros(numVerts, 2);
vertexCoordsFit = zeros(numVerts, 2);

% Calculate the 'uncalibrated' vertex coordinates
for y=1:yLoomSize
    for x=1:xLoomSize
        index = ((y-1) * xLoomSize) + x;
        xCoord = (x-1) * xStep;
        yCoord = (y-1) * yStep;
        vertexCoords(index, :) = [xCoord yCoord];
        %fprintf('(%f, %f)\n', vertexCoords(index, 1), vertexCoords(index, 2));
    end
end

% Some debug plots, if requested:
if showCalibOutput
    figure(9)
    hold off
    plot(scal.XCALIBDOTS, scal.YCALIBDOTS, 'b.')
    hold on
    plot(scal.SELECTXCALIBDOTS, scal.SELECTYCALIBDOTS, 'r.')
    size(scal.XCALIBDOTS)
    size(scal.YCALIBDOTS)
    size(scal.SELECTXCALIBDOTS)
    size(vertexCoords(:,1))
    size(vertexCoords(:,2))
end

% Fit mesh vertices to locations of calibrated points - Use Matlabs
% griddata fitting and interpolation routine:
if ~IsOctave
    % Matlab: Use 'v4' method - the interpolation method used by Matlab V4
    % for interpolation:
    vertexCoordsFit(:,1)= griddata(scal.XCALIBDOTS, scal.YCALIBDOTS, scal.SELECTXCALIBDOTS, vertexCoords(:,1), vertexCoords(:,2), 'v4'); %#ok<*GRIDD>
    vertexCoordsFit(:,2)= griddata(scal.XCALIBDOTS, scal.YCALIBDOTS, scal.SELECTYCALIBDOTS, vertexCoords(:,1), vertexCoords(:,2), 'v4');
else
    % Octave: griddata() is also supported by GNU/Octave, but the 'v4' method is
    % only supported by Matlab, not by Octave. 'Cubic' isn't supported
    % either, so the best we can do is to use the default 'linear' method.
    % This is problematic as it creates different results when running on
    % Octave vs. Matlab:
    % Ok - Actually it doesn't work at all on Octave, because the relevant
    % implementation of Octave's griddata() seems to be quite buggy :-(
    % TODO FIXME: Is this still true? Octave 3.4 supports 'cubic' and maybe
    % 'linear' has been fixed already?
    vertexCoordsFit(:,1)= griddata(scal.XCALIBDOTS, scal.YCALIBDOTS, scal.SELECTXCALIBDOTS, vertexCoords(:,1), vertexCoords(:,2));
    vertexCoordsFit(:,2)= griddata(scal.XCALIBDOTS, scal.YCALIBDOTS, scal.SELECTYCALIBDOTS, vertexCoords(:,1), vertexCoords(:,2));
end

% Some debug plots, if requested:
if showCalibOutput
    figure(10)
    hold off
    plot(vertexCoords(:,1), vertexCoords(:,2), 'b.')
    hold on
    plot(vertexCoordsFit(:,1), vertexCoordsFit(:,2), 'r.')

    figure(11)
    hold off
    plot(scal.SELECTXCALIBDOTS, scal.SELECTYCALIBDOTS, 'o')
    hold on
    plot(vertexCoordsFit(:,1), vertexCoordsFit(:,2), 'r.')
end

% Fit the coordinates to the calibrated space (values are in pixels)
%    [vertexCoordsFit(:,1) vertexCoordsFit(:,2)] = BFbvlFitCoords(vertexCoords(:,1), ...
%        vertexCoords(:,2), xFitCoef_R, yFitCoef_R);

% Compute final vertex- and texcoords. Need to swap y-positions upside-down
% as our internal vertex/texcoord assignment is upside down wrt. original
% Banks lab calibration:
vertexCoords(:,2) = windowHeight - vertexCoords(:,2);
[xyzcalibpos, xytexcoords] = BVLGeneratetextcoord(yLoomSize, xLoomSize, vertexCoords, vertexCoordsFit, showCalibOutput);

% Done. Return results:
return;


function [xyzcalibpos, xytexcoords]=BVLGeneratetextcoord(yLoomSize, xLoomSize, vertexCoords, vertexCoordsFit, showCalibOutput)
% Internal helper routine: Called by BVLComputeWarpMesh.m.
% Implements the geometric display calibration and undistortion procedure
% developed by the Banks Vision Lab at UC Berkeley.
%
% History:
% 02/17/08  Derived from BFGeneratetextcoord.m with minimal modifications. (MK)
%

numVerts = (xLoomSize-1) * (yLoomSize-1) * 4;

xyzcalibpos = zeros(1, numVerts*2);
xytexcoords = zeros(1, numVerts*2);

xtemp = zeros(1, numVerts);
ytemp = zeros(1, numVerts);

xverts= zeros(1, numVerts);
yverts= zeros(1, numVerts);

vectaddress=0;

if showCalibOutput
    figure(100);
    axis ij;
    hold on;
    plot(vertexCoords(:,1), vertexCoords(:,2), 'r.');
    plot(vertexCoordsFit(:,1), vertexCoordsFit(:,2), 'b.');
end

for y=1:(yLoomSize-1)
    for x=1:(xLoomSize-1)
        index = ((y-1) * xLoomSize) + x;

        vectaddress=vectaddress+1;

        xtemp(vectaddress)  = vertexCoords(index, 1);  %LowerLeftXTex
        ytemp(vectaddress)  = vertexCoords(index, 2);  %LowerLeftYTex

        xverts(vectaddress) = vertexCoordsFit(index, 1);  % Lower left fit coord  x
        yverts(vectaddress) = vertexCoordsFit(index, 2);  % lower left fit coord  y

        vectaddress=vectaddress+1;

        xtemp(vectaddress)  = vertexCoords(index+1, 1);  %LowerrightXTex
        ytemp(vectaddress)  = vertexCoords(index+1, 2); %LowerrightYTex

        xverts(vectaddress) = vertexCoordsFit(index+1, 1);   %lower right fit coord  x
        yverts(vectaddress) = vertexCoordsFit(index+1, 2);   %lower right fit coord  y

        vectaddress=vectaddress+1;

        xtemp(vectaddress)  = vertexCoords(index+xLoomSize+1, 1);  %UpperRightXTex
        ytemp(vectaddress)  = vertexCoords(index+xLoomSize+1, 2);  %UpperRightYTex

        xverts(vectaddress) = vertexCoordsFit(index+xLoomSize+1, 1);  %Upper right fit coord  x
        yverts(vectaddress) = vertexCoordsFit(index+xLoomSize+1, 2);  %upper right fit coord  y

        vectaddress=vectaddress+1;

        xtemp(vectaddress)  = vertexCoords(index+xLoomSize, 1);  %UpperLeft X Tex
        ytemp(vectaddress)  = vertexCoords(index+xLoomSize, 2);  %Upperleft Y tex

        xverts(vectaddress) = vertexCoordsFit(index+xLoomSize, 1);   %upper left fit coord   x
        yverts(vectaddress) = vertexCoordsFit(index+xLoomSize, 2);   %upper left fit coord   y
    end
end

xyzcalibpos(1:2:end) = xverts;
xyzcalibpos(2:2:end) = yverts;
xytexcoords(1:2:end) = xtemp;
xytexcoords(2:2:end) = ytemp;

return;

% --- End of Helper routines for setup of the calibration method 'BVLDisplayList' ---

% --- Helper routines for setup of the calibration method 'CSVDisplayList' ---

function [xyzcalibpos, xytexcoords] = CSVComputeWarpMesh(windowWidth, windowHeight, scal, showCalibOutput)
% [xyzcalibpos, xytexcoords] = CSVComputeWarpMesh(windowWidth, windowHeight, scal, showCalibOutput)
%
% Use the calibration information stored in 'scal', together with the
% current 'windowWidth' and 'windowHeight' of the onscreen window to
% generate vectors of vertex coordinates and texture coordinates for a mesh
% that performs the proper "display undistortion".
%
% History:
% 07/26/12  Derived from BVLComputeWarpMesh.m with major simplifications. (MK)
%

% Generate the calibration vertices here so we only do this once.
xLoomSize = size(scal.vcoords, 2);
yLoomSize = size(scal.vcoords, 1);
numVerts = xLoomSize * yLoomSize;

% vertexCoords = Nx2 array, N rows of [x y] pairs. Row-Major format encoding.
textureCoords = zeros(numVerts, 2);
vertexCoords  = zeros(numVerts, 2);

% Parse the matrices passed in scal and rearrange them to the format of the
% vertexCoords and textureCoords vectors: Scanning is row-major order.
% We also scale all positions with window width and height, as the scal
% matrices contain normalized coordinates in 0.0 - 1.0 range for display
% width/height of a "unit display". vcoords can exceed that range or be
% negative - they are assigned to positions outside the framebuffer.
for y=1:yLoomSize
    for x=1:xLoomSize
        index = ((y-1) * xLoomSize) + x;
        vertexCoords(index, :)  = [scal.vcoords(y, x, 1) * windowWidth, scal.vcoords(y, x, 2) * windowHeight];
        textureCoords(index, :) = [scal.tcoords(y, x, 1) * windowWidth, scal.tcoords(y, x, 2) * windowHeight];
    end
end

% Compute final vertex- and texcoords. Need to swap y-positions upside-down
% as our internal vertex/texcoord assignment is upside down wrt. original
% calibration:
textureCoords(:,2) = windowHeight - textureCoords(:,2);

% Some debug plots, if requested:
if showCalibOutput
    figure;
    hold on;
    axis ij;
    plot(textureCoords(:,1), textureCoords(:,2), 'b.');
    plot(vertexCoords(:,1), vertexCoords(:,2), 'r.');
    hold off;
end

% textureCoords are regularly spaced texture 2D coordinates.
% vertexCoords are irregularly placed vertex 2D coordinates.
[xyzcalibpos, xytexcoords] = BVLGeneratetextcoord(yLoomSize, xLoomSize, textureCoords, vertexCoords, showCalibOutput);

% Done. Return results:
return;

% --- End of helper routines for setup of the calibration method 'CSVDisplayList' ---
