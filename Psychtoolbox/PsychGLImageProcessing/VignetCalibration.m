function s = VignetCalibration(filename, screenId)
% s = VignetCalibration([filename] [, screenId])
%
% Vignetted Luminance calibration procedure for undistortion of distorted
% display luminance.
%
% This code was created by Jorrit S. Montijn of the Laboratoire Psychologie
% de la Perception at Universite Paris Descartes.
%
% This function allows the user to create a calibration for a vignetted
% luminance distortion. It enables the user to choose from either a
% Gaussian or Exponential approximation of the luminance vignetting. Users
% can change the position and steepness of the luminance fall-off in
% horizontal and vertical direction independently. It is also possible to
% adjust the minimal luminance.
%
% Once the calibration procedure is finished, the gain matrix (and other
% possibly useful parameters and variables) are saved to disk for use by a
% vignetting compensation program (see 'VignettingCorrectionDemo' for an
% example of how to do this).
% 
% How to use:
% 
% Start the script, providing a filename and screennumber:
%
% 'filename' is the name of the file to which calibration results will be 
% saved. If none is provided, a default name will be used. The default file
% will be stored in the directory returned by PsychtoolboxConfigDir('ShadingCalibration').
%
% 'screenId' is the optional index of the screen that you want to
% calibrate. If only one screen is present, it is not necessary to include
% this parameter.
%
% When started, the program will show either a gaussian or exponential
% approximation of the required gain function. You can adjust the
% parameters so that the screen will look approximately isoluminant. It is
% recommended to use a photometer, since human eyes are far from perfect when
% doing this.
%
% You can use the following keys to issue commands:
%
% Escape	-	use this key to save and exit
% Tab		-	show/hide help text
% Space		-	change the parameter you wish to adjust
% '<' and '>'	-	use these keys to change the increment size of the
%					parameter adjustments
% Arrow Keys	-	adjust parameters (depending on selected parameter)
%				1. Luminance/Function type
%				2. Width/Steepness of function
%				3. Position of centre
%

% History:
%
% 16.03.2010     Interactive Vignetting Calibration function created
%			     by Jorrit S. Montijn.
%
% 18.03.2010  mk Minimal cleanup and merge into distribution.
% 27.03.2010  mk Switch to multiplicative inverse. Use gainCorrection of a
%                50% grey background via PsychColorCorretion(). This allows
%                for > 1.0 gain values.
% 19.04.2010  jm Removed Lmin restrictions allowing for gain values to take
%				 any value (-inf < Lmin < inf), added display of 
%				 current parameter values in on-screen help text
%

	% Check for proper installation:
	PsychDefaultSetup(1);

	%create variable if needed
    if ~exist('screenId','var')
        screenId = [];
    end

    % Fetch connected screens:
    screens = Screen('Screens');

    %if no screen nr was provided
    if isempty(screenId)
        %ask which one to calibrate if more than 1 is present
        if length(screens)>1,
            screenId = str2double(input(['Which screen should we calibrate? [' num2str(screens) '] '],'s'));
        else
            screenId = screens;
        end
    end

    if ~ismember(screenId, screens)
        error('No such screen. Please restart program with valid screenid.');
    end
    
	%put screenId number into structure
	s.scrNr = screenId;

	%fetch screen resolution
	[w, h] = Screen('WindowSize', screenId);

	%create variable if needed
	if ~exist('filename', 'var')
		filename = [];
	end

	%if no filename was provided
	if isempty(filename)
		%assign filename
		filename = [PsychtoolboxConfigDir('ShadingCalibration') 'VignetCalibration' sprintf('_%i_%i_%i', screenId, w, h) '.mat'];
		fprintf('\nNo filename defined... Using default\n');
	end

	%show name of output file and check to overwrite
	fprintf('Name of outputfile: %s\n\n', filename);
	if exist(filename, 'file')
		answer = input('File already exists. Overwrite? [y/n]','s');
		if ~strcmp(lower(answer), 'y') %#ok<STCI>
			fprintf('\n\nCalibration aborted. Please restart program and specify another filename.\n\n');
			return;
		end
	end

	%disable synctests
	oldsynclevel = Screen('Preference', 'SkipSyncTests', 2);
	%change verbosity
	oldverbosity = Screen('Preference', 'Verbosity', 1);

	% Open window on desired screenId:
    % Use imaging pipeline for visual processing and display:
    PsychImaging('PrepareConfiguration');

    % Request per-pixel 2D gain correction for display:
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'GainMatrix');

    % Open window with gray background (128 aka 50% gray), return handle 'win':
    [s.window, rect] = PsychImaging('OpenWindow', s.scrNr, 128);
    
	%change text appearance
	Screen('TextSize',s.window, 18);
	Screen('TextStyle', s.window, 1+2);
	Screen('TextColor', s.window, [255 0 0]);

	%put rect into structure
	s.rect = rect;

	try
		%hide cursor and disable input
		HideCursor;
		ListenChar(2);

		%initialize parameters
		[beta,xGrid,yGrid] = InitVignetCalibration(s.rect);

		%put values into structure
		s.quitNow = 0;
		s.beta = beta;
		s.xGrid = xGrid;
		s.yGrid = yGrid;

		%start calibration
		s = CalibrateVignetNow(s);

		%collect variables from structure for easy access
		gainMatrix = s.textureGrid; %#ok<NASGU>
		parameters = s.beta; %#ok<NASGU>
		functionType = s.functionType; %#ok<NASGU>
		
        % Strip redundant data from 's':
        bs = s;
        s.xGrid = [];
        s.yGrid = [];
        s.textureGrid = [];
        
		%save variables and structure to file
		save(filename, 'gainMatrix', 'parameters', 'functionType', 's', '-mat', '-V6');

        % Restore unstripped 's' before return:
        s = bs;
        
		fprintf('The calibration is finished! \n\n');
		fprintf('Your calibration matrix, parameters and function type have been saved to\n%s\n', filename);
		fprintf('Also included in the file is a structure that contains most settings you used.\n');
		fprintf('For an example on how to use the generated gain matrix to reduce vignetting,\n');
		fprintf('please take a look at the PsychToolbox demo ''VignettingCorrectionDemo''\n\n');
	catch
		%show cursor
		ShowCursor('Arrow');

		%reset screen values
		Screen('CloseAll');
		Screen('Preference', 'SkipSyncTests', oldsynclevel);
		Screen('Preference', 'Verbosity', oldverbosity);

		%reenable listenchar
		ListenChar(0);

		%rethrow error
		psychrethrow(psychlasterror);
	end

	%reenable listenchar
	ListenChar(0);

	%restore cursor
	ShowCursor('Arrow');

	%reset screen values
	Screen('CloseAll');
	Screen('Preference', 'SkipSyncTests', oldsynclevel);
	Screen('Preference', 'Verbosity', oldverbosity);
	
	return;
end

function [beta,xGrid,yGrid] = InitVignetCalibration(rect)

	rLeft = rect(1);
	rTop = rect(2);
	rRight = rect(3);
	rBot = rect(4);

	xPixelSpace = rLeft+1:rRight; %horizontal pixel space
	yPixelSpace = rTop+1:rBot; %vertical pixel space

	xPixNum = length(xPixelSpace); %number of pixels on x-axis
	yPixNum = length(yPixelSpace); %number of pixels on y-axis

	Lmin = 0; %minimal luminance level
	wX = xPixNum/2; %initial horizontal gaussian width is quarter of screen
	wY = yPixNum/2; %initial vertical gaussian width is quarter of screen
	muX = xPixNum/2; %initial horizontal gaussian position is at middle of screen
	muY = yPixNum/2; %initial vertical gaussian position is at middle of screen

	%put parameters into beta variable
	beta(1) = Lmin;
	beta(2) = wX;
	beta(3) = wY;
	beta(4) = muX;
	beta(5) = muY;

	%create meshgrids from pixel spaces
	[xGrid,yGrid] = meshgrid(xPixelSpace,yPixelSpace);
end


function s = CalibrateVignetNow(s)

	%define initial variables
	useFunction = 1;
	drawText = 1;
	paramNow = 1;
	beta = s.beta;
	xGrid = s.xGrid;
	yGrid = s.yGrid;
	switchVar = [2 1];
	switchVar3 = [2 3 1];
	redrawTexture = 1;
	functionList = {'Gaussian','Exponential'};
	paramList = {'Minimum Luminance','Width','Position'};

	%create Lmin variables
	Lmin = beta(1);
	LminIncrement = 0.1;

	%create width variables
	wX = beta(2);
	wY = beta(3);
	wIncrement = wX/20;

	%create position variables
	muX = beta(4);
	muY = beta(5);
	muIncrement = muX/20;

	while s.quitNow == 0
		%if texture whould be (re)drawn
		if redrawTexture == 1
			%reset redraw var
			redrawTexture = 0;

			%create grid from parameters
			beta(1) = Lmin;
			beta(2) = wX;
			beta(3) = wY;
			beta(4) = muX;
			beta(5) = muY;
			if useFunction == 1
				%use gaussian approximation
				textureGrid = gaussianGrid(beta, xGrid, yGrid);
			else
				%use exponential approximation
				textureGrid = exponentialGrid(beta,xGrid,yGrid);
			end

            % Apply gain correction matrix to display: This will take effect at the
            % next Screen('Flip'):
            PsychColorCorrection('SetGainMatrix', s.window, textureGrid);
		end

		%if text should be drawn
		if drawText == 1
			%draw help
			s.functionType = functionList{useFunction};
			if paramNow == 1
				thisIncrement = num2str(LminIncrement);
				upDown = 'increase/decrease min Luminance factor';
				curValUD = num2str(beta(1));
				leftRight = 'Shift to other function';
				curValLR = s.functionType;
			elseif paramNow == 2
				thisIncrement = num2str(wIncrement);
				upDown = 'increase/decrease vertical sigma';
				curValUD = num2str(beta(3));
				leftRight = 'decrease/increase horizontal sigma';
				curValLR = num2str(beta(2));
			else
				thisIncrement = num2str(muIncrement);
				upDown = 'move centre up/down';
				curValUD = num2str(beta(5));
				leftRight = 'move centre left/right';
				curValLR = num2str(beta(4));
			end
			
			text = ['You can use the following commands:\n\n'...
				'TAB: show/hide this text\n'...
				'Escape: save&exit\n'...
				'Space: switch parameters :: Current: ' paramList{paramNow} '\n'...
				'"<" and ">": decrease/increase increment: ' thisIncrement '\n'...
				'Up/Down: ' upDown ' :: Current: ' curValUD '\n'...
				'Left/Right: ' leftRight ' :: Current: ' curValLR '\n'];
			DrawFormattedText(s.window,text,20,20);
		end

		%flip screen
		Screen('Flip',s.window);

		%check for input
		inputReceived = 0;
		keyPressed = 0;
        while inputReceived == 0 && s.quitNow == 0
            KbWait; %wait till keypress
            [keyIsDown,secs,keyCode]=KbCheck;
            if keyIsDown
                if keyPressed == 0
                    keyPressed = 1;
                    inputReceived = 1;
                    if keyCode(KbName('Escape'))
                        %done! Put output texture  and parameters in
                        %structure to return to master function and break
                        %loop

                        s.textureGrid = textureGrid;
                        s.beta = beta;
                        s.quitNow = 1;
                        break;
                    elseif keyCode(KbName('space'))
                        %shift adjustment parameters
                        %1: strength (Lmin)
                        %2: width (wX and xY)
                        %3: position (muX and muY)
                        paramNow = switchVar3(paramNow);
                    elseif keyCode(KbName('tab'))
                        %toggle show/hide text
                        drawText = switchVar(drawText);
                    else
                        if paramNow == 1
                            %Lmin
                            valueX = Lmin;
                            valueY = Lmin;
                            increment = LminIncrement;
                        elseif paramNow == 2
                            %width
                            valueX = wX;
                            valueY = wY;
                            increment = wIncrement;
                        elseif paramNow == 3
                            %position
                            valueX = muX;
                            valueY = muY;
                            increment = muIncrement;
                        end

                        if keyCode(KbName('LeftArrow'))
                            %adjust left
                            redrawTexture = 1;
                            valueX = valueX - increment;
                            if paramNow == 1
                                useFunction = switchVar(useFunction);
                            end
                        elseif keyCode(KbName('RightArrow'))
                            %adjust right
                            redrawTexture = 1;
                            valueX = valueX + increment;
                            if paramNow == 1
                                useFunction = switchVar(useFunction);
                            end
                        elseif keyCode(KbName('UpArrow'))
                            %adjust up
                            redrawTexture = 1;
                            if paramNow == 3
                                valueY = valueY - increment;
                            else
                                valueY = valueY + increment;
                            end
                        elseif keyCode(KbName('DownArrow'))
                            %adjust down
                            redrawTexture = 1;
                            if paramNow == 3
                                valueY = valueY + increment;
                            else
                                valueY = valueY - increment;
                            end
                        elseif keyCode(KbName(',<'))
                            %decrease increment
                            increment = increment / 2;
                        elseif keyCode(KbName('.>'))
                            %increase increment
                            increment = increment * 2;

                        else
                            %false alarm: no input received, just a random keypress
                            inputReceived = 0;
                        end

                        if paramNow == 1
                            %Lmin
                            LminIncrement = increment;
                            Lmin = valueY;
                        elseif paramNow == 2
                            %width
                            wIncrement = increment;
                            wX = valueX;
                            wY = valueY;
                        elseif paramNow == 3
                            %position
                            muIncrement = increment;
                            muX = valueX;
                            muY = valueY;
                        end
                    end
                end
            else
                keyPressed = 0;
            end
        end
        KbReleaseWait;
	end
end


function gaussGrid = gaussianGrid(beta, xGrid, yGrid)
	
	Lmin = beta(1); %minimal luminance level (0-1)
	wX = beta(2); %width on x-axis (pixels)
	wY = beta(3); %width on y-axis (pixels)
	muX = beta(4); %centre on x-axis (pixels)
	muY = beta(5); %centre on y-axis (pixels)

	gaussX = exp(- (((xGrid-muX).^2) / ((2*wX).^2) ) ); %create horizontal gaussian
	gaussY = exp(- (((yGrid-muY).^2) / ((2*wY).^2) ) ); %create vertical gaussian

	comb = gaussX.*gaussY; %combine into one grid

    % Invert grid: Multiplicative inverse should compensate for
    % multiplicative attenuation of display. Add offset Lmin to shift the
    % gainfield around: The final matrix stores floating point gain values
    % greater or equal to 0.0:
	gaussGrid = (1 ./ comb) + Lmin;

    fprintf('Mingain = %f, Maxgain = %f\n', min(min(gaussGrid)), max(max(gaussGrid)));
end


function expoGrid = exponentialGrid(beta, xGrid,yGrid)

	Lmin = beta(1); %minimal luminance level (0-1)
	xF = beta(2); %width on x-axis (pixels)
	yF = beta(3); %width on y-axis (pixels)
	muX = beta(4); %centre on x-axis (pixels)
	muY = beta(5); %centre on y-axis (pixels)


	expX = min(exp(xGrid/xF-muX/xF),exp(-xGrid/xF+muX/xF) ); %horizontal part
	expY = min(exp(yGrid/yF-muY/yF),exp(-yGrid/yF+muY/yF) ); %vertical part

	comb = expX.*expY; %combine into one grid

    % Invert grid: Multiplicative inverse should compensate for
    % multiplicative attenuation of display. Add offset Lmin to shift the
    % gainfield around: The final matrix stores floating point gain values
    % greater or equal to 0.0:
	expoGrid = (1 ./ comb) + Lmin;
	
    fprintf('Mingain = %f, Maxgain = %f\n', min(min(expoGrid)), max(max(expoGrid)));
end
