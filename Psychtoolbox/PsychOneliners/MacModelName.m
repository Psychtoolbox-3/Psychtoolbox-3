function modelName=MacModelName
% modelName=MacModelName
% Return the model name of the Macintosh. 
%
% OSX: Get the computer identifier as reported by open firmware.  In the
% OSX Psychtoolbox that is provided in the struct returned by
% Screen('Computer'). Convert that identifier to Apple's marketing
% name for the computer.  To do that we use the mapping of identifier to
% product names found in:
% /System/Library/SystemProfiler/SPPlatformReporter.spreporter/Contents/
%  Resources/English.lproj/Localizable.strings.  
% There is no official Apple-recommended way to uncover the marketing
% name. However, this method works and should be reasonably robust.
% 
% See also: Gestalt, Screen('Computer?'), OSName, AppleVersion

% HISTORY
%
% 2004    awi Wrote it.
% 1/29/05 dgp Cosmetic.
% 3/14/05 dgp Fixed handling of model names that end in "\"".
% 3/5/06  awi Added note explaining new bug and how to repair.
% 5/29/06  mk Added dumb fix for endian issue on new Intel-Macs.
% 9/20/09  mk Don't fail with name 'Unknown' but output raw model name then. 

% NOTES
%
% As of 10.4.5 Apple has changed the file format in which System Profiler
% stores the map of system firmware IDs to computer model names.  The map
% is now a property list stored in this file:
% /System/Library/SystemProfiler/SPPlatformReporter.spreporter/Contents/Resources/SPMachineTypes.plist
%
% The easiest and most structured way to adapt the PTB to this change would
% be to provide a mex function which reads a property list file and returns
% a MATLAB structure holding the contents of the property list. Because
% functions provided by MacOS would do the work of parsing the XML property
% list file, that should be a simple program to write.  Otherwise, there
% might even already exist XML parsers for MATLAB.  


if IsOSX
    % The mapping file is in Unicode format, which MATLAB does not understand.  We can't treat
    % the content as characters without doing some work first.  It would be
    % nice if we could use built-in MATLAB routines such as fgetl, but no dice.
    persistent modelNameCache
    
    if ~isempty(modelNameCache)
        modelName=modelNameCache;
        return
    end

    mappingFileName='/System/Library/SystemProfiler/SPPlatformReporter.spreporter/Contents/Resources/English.lproj/Localizable.strings';
    fp=fopen(mappingFileName,'r');
    if fp==-1
        error(['Failed to open name mapping file: ' mappingFileName]);
    end
    fileWords=fread(fp,inf,'uint16');
    fclose(fp);
    % Unicode files start with hex feff or ffef to indicate byte order.  
    % Confirm that it's feff, and then strip it off.
    if fileWords(1)~=hex2dec('feff');
        % Convert fileWords by first casting it to double type, then
        % kind of endian-swapping. We throw away the least significant
        % byte value and store the most significant byte as value. This
        % would discard any former low-bytes. Dont know if this will work.
        fileWords=floor(double(fileWords)/256);
    end
    fileWords=fileWords(2:end);
    fileCharsRaw=char(fileWords);
    % Strip out comments.
    commentsOnFlag=0;
    fileChars='';
    i=1;
    while i <= length(fileCharsRaw)-2 %we don't care about comment toggles in the last two characters.
        if commentsOnFlag==0 && fileCharsRaw(i) == '/' && fileCharsRaw(i+1) == '*'
            commentsOnFlag=1;
            i=i+2;
        elseif commentsOnFlag==1 && fileCharsRaw(i) == '*' && fileCharsRaw(i+1) == '/'
            commentsOnFlag=0;
            i=i+2;
        else
            if ~commentsOnFlag
                fileChars(end+1)=fileCharsRaw(i);
            end
            i=i+1;
        end
            
    end
    % Break the file into lines.
    lineBreakChar=char(10); 
    currentLineIndex=1;
    currentCharIndex=1;
    fileLines{currentLineIndex}='';
    while currentCharIndex <= length(fileChars)
        if fileChars(currentCharIndex)==lineBreakChar;
            currentLineIndex=currentLineIndex+1;
            fileLines{currentLineIndex}='';
        else
            fileLines{currentLineIndex}(end+1)=fileChars(currentCharIndex);
        end
        currentCharIndex=currentCharIndex+1;
    end
    % Strip out blank lines.
    blanckLineIndices=find(streq(fileLines, ''));
    contentLineIndices=setdiff(1:length(fileLines), blanckLineIndices);
    contentLines={fileLines{contentLineIndices}};
    % We now have a cell array of content lines.  Consider the stuff we did up to now to be preprocessing.
    % Ad hoc parsing begins here.
    %
    % Use lines left and right of the equals sign to create a dictionary.
    KEY=1;
    ENTRY=2;
    acc=1;
    dictV1={};
    for i=1:length(contentLines)
        equalLocations=find(contentLines{i}=='=');
        if length(equalLocations) == 1
            dictV1{acc,KEY}=contentLines{i}(1:equalLocations-1);
            dictV1{acc,ENTRY}=contentLines{i}(equalLocations+1:end);
            acc=acc+1;
        elseif length(equalLocations) > 1
            error('MacModelName encountered an entry with more than one equals sign and parsing failed');
        % elseif  length(equalLocations) < 1
            % should not happen but ignore it.    
        end
    end
    % Intereresting lines look like this:
    %   "PowerBook5,4" = "PowerBook G4 15\"";
    % Uninteresting lines are distinguished by the absence of quotes on the
    % left:
    %   number_processors = "Number Of CPUs";
    % Filter out uninteresting lines
    %
    dictV2={};
    acc=1;
    for i=1:length(dictV1);
        if find(dictV1{i,KEY}=='"') > 0;
            dictV2(acc,:)=dictV1(i,:);
            acc=acc+1;
        end
    end
    % Get the identifier and look it up in the dictionary.  Note that the
    % dictionary key contains more characters than the platform identifier,
    % such as spaces and quotes, so we check only for membership, not
    % identity.
    c=Screen('Computer');
    foundEntry=0;
    for i=1:length(dictV2);
        if strfind(dictV2{i,KEY}, c.hw.model)
            foundEntry=1;
            rawName=dictV2{i,ENTRY};
            break;
        end
    end
    if foundEntry
        surroundingQuoteIndices=strfind(rawName,'"');
        modelName=rawName(surroundingQuoteIndices(1)+1:surroundingQuoteIndices(end)-1);
    else
        modelName=c.hw.model;
    end
    if modelName(end-1)=='\'
        modelName=modelName([1:end-2 end]);
    end
    modelNameCache=modelName;
end
