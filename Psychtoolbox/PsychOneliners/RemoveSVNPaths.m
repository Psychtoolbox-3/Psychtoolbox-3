function newPathList = RemoveSVNPaths(pathList)
% newPathList = RemoveSVNPaths(pathList)
% Removes any .svn paths from the pathList.  If no pathList is specified,
% then the program sets pathList to the result of the 'path' command.  This
% function returns a 'pathsep' delimited list of paths omitting the .svn
% paths.

% History:
% 14.07.06 Written by Christopher Broussard.
% 25.07.06 Modified to work on M$-Windows and GNU/Octave as well (MK).
% 31.05.09 Adapted to fully work on Octave-3 (MK).

% If no pathList was passed to the function we'll just grab the one from
% Matlab.
if nargin ~= 1
    % Grab the path list.
    pathList = path;
end

try
    % We do the .svn path removal in a try-catch block, because some of the
    % functions used inside this block are not available in Matlab-5 and
    % GNU/Octave. Our catch - block provides fail-safe behaviour for that
    % case.
    
    % Break the path list into individual path elements.
    if IsOctave
        pathElements = strsplit(pathList,pathsep,true);
    else
        pathElements = textscan(pathList, '%s', 'delimiter', pathsep);
        pathElements = pathElements{1}.';
    end
    
    % Look at each element from the path.  If it doesn't contain a .svn folder
    % then we add it to the end of our new path list.
    qNotSVN = cellfun(@isempty,strfind(pathElements,[filesep '.svn']));
    pathElements = pathElements(qNotSVN);
    
    if ~isempty(pathElements)
        % generate new pathList
        pathElements = [pathElements; repmat({pathsep},1,length(pathElements))];
        newPathList  = [pathElements{:}];
        
        % drop last separator
        newPathList(end) = [];
    end
catch
    % Fallback behaviour: We fail-safe by simply returning the unmodified
    % pathList. No .svn paths removed, but the whole beast is still
    % functional.
    newPathList = pathList;
end
