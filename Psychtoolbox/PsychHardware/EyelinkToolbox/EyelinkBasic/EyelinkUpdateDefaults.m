function EyelinkUpdateDefaults(el)
% USAGE = EyelinkUpdateDefaults(el) 
%
% This function passes changes to the calibration defaults structure
% to the callback function. To be called after EyelinkInitDefaults if any
% changes are made to the output structure of that function, e.g.
%    el=EyelinkInitDefaults(window);
%    el.backgroundColour = BlackIndex(window);
%    EyelinkUpdateDefaults(el);
%
% 27-1-2011 NJ created 
% 19-12-2012 IA Fix hardcoded callback 

if ~isempty(el.callback) && exist(el.callback,'file')
    if ~isempty(el.calImageTargetFilename)
        if exist(el.calImageTargetFilename, 'file')
            el.calImageInfo = imfinfo(el.calImageTargetFilename); % Get image file info
            el.calImageData = imread(el.calImageTargetFilename); % Read image from file
            el.calImageTexture = Screen('MakeTexture', el.window, el.calImageData); % Convert image file to texture
            if ~strcmpi(el.calTargetType, 'image') 
                warning(sprintf([ ...
                    'EyelinkToolbox -- ''el.calImageTargetFilename'' is configured but will not be used until\n' ...
                    '\tel.calTargetType = ''image''\n' ...
                    'is also set and \n' ...
                    '\tEyelinkUpdateDefaults(el)\n' ...
                    'is called again thereafter.\n\n']));
            end
        else
            warning(sprintf([ ...
                'EyelinkToolbox -- File Not Found:\n', ...
                '\tel.calImageTargetFilename = %s\n\n'], ... 
                el.calImageTargetFilename));
        end
    end

    if ~isempty(el.calAnimationTargetFilename)
        if exist(el.calAnimationTargetFilename, 'file')
            if strcmpi(el.calTargetType, 'video')
                if el.targetbeep ~= 0
                   warning(sprintf([ 'EyelinkToolbox - ''el.calAnimationTargetFilename'' set for video, but\n' ...
                       'el.targetbeep not set == 0 and may cause playback issues (freezing).\n']));
                end
                
                if el.feedbackbeep ~= 0
                   warning(sprintf([ 'EyelinkToolbox - ''el.calAnimationTargetFilename'' set for video, but\n' ...
                       'el.feedbackbeep not set == 0 and may cause playback issues (freezing).\n']));
                end
            else
                warning(sprintf([ ...
                    'EyelinkToolbox -- ''el.calAnimationTargetFilename'' is configured but will not be used until\n' ...
                    '\tel.calTargetType = ''video''\nis also set and \n\tEyelinkUpdateDefaults(el)\nis called again thereafter.\n\n']));
            end
        else
            warning(sprintf([ ...
                'EyelinkToolbox -- File Not Found:\n', ...
                '\tel.calAnimationTargetFilename = %s\n\n'], ...
                el.calAnimationTargetFilename));
        end
    end

    % evaluate the callback function with the new el defaults
    feval(el.callback, el); 
end
end
