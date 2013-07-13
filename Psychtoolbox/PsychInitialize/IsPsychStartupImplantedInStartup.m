function isImplanted = IsPsychStartupImplantedInStartup
% isImplanted = IsPsychStartupImplantedInStartup
%
% Return TRUE if both the startup.m file is on the MATLAB path and if it
% includes a call to PsychStartup.
%
% If you have your own startup.m file which you want to use, then place it
% ahead of the Psychtoolbox startup.m file on MATLAB's search path and add
% to your file the line:
%   PsychStartup; % Setup Psychtoolbox.
% _________________________________________________________________________


commandString='PsychStartup';
commentChar='%';

pth=which('startup.m');
if isempty(pth)
    isImplanted=0;
else
    %check that the startup file contains the PsychStartup call and don't be
    %fooled by commented-out linese. we implement these as a series of filters.
    
    %read the file contents and breack into lines
    f=fopen(pth);
    cDoubles=fread(f);
    fclose(f);
    cChars=char(cDoubles');
    cLines=BreakLines(cChars);
    %find the indices of those lines which contain the the string
    %"PsychStartup".  These are all candidates but some
    % might just be comments.
    occuranceIndices=[];
    for i =1:length(cLines)
        if ~isempty(strfind(cLines{i}, commandString))
            occuranceIndices=[occuranceIndices i]; %#ok<*AGROW>
        end
    end
    %build a cell array of the line fragments before the comment character.
    nonCommentFragments={};
    if length(occuranceIndices) < 1
        isImplanted = 0;
        return;
    end
    
    for i=length(occuranceIndices)
        commentCharIndices=find(cLines{occuranceIndices(i)}==commentChar);
        if isempty(commentCharIndices);
            %add the line to the list of candidate fragments if there is no
            %comment character
            nonCommentFragments{i}=cLines{occuranceIndices(i)};
        else
            %add to the list of candidate fragments the portion of the line
            %ahead of the comment character if there is one.
            temp=cLines{occuranceIndices(i)};
            nonCommentFragments{i}=temp(1:commentCharIndices(1));
        end
    end
    %search the list of candidate fragments for the command string.  This
    %is imperfect because it can be fooled if the command is a substring of
    %a longer command.  But its probably good enough.
    numFoundCalls=0;
    foundCallsLineIndices=[];
    for i=length(occuranceIndices)
        foundIndices=strfind(nonCommentFragments{i},commandString );
        if ~isempty(foundIndices)
            numFoundCalls = numFoundCalls+1;
            foundCallsLineIndices = [foundCallsLineIndices i];
        end
    end
    if numFoundCalls > 1
        % we can't print multiline warnings using warning so just print the
        % the word "warning"
        fprintf('Warning: The startup.m script:\n');
        fprintf(['\t' pth '\n']);
        fprintf(['seems to call "' commandString '" ' int2str(numFoundCalls) ' times.  Only one call is necessary\n']);
    end
    isImplanted = numFoundCalls;
end
