function newPathList = RemoveMatchingPaths(pathList, matchString)
% newPathList = RemoveMatchingPaths(pathList, matchString)
%
% Removes any paths that contain the given matchString from the pathList.
% If no pathList is specified, then the program sets pathList to the result
% of the 'path' command.  This function returns a 'pathsep' delimited list
% of paths, omitting the paths that contained the given matchString.

% History:
% 30.05.13 Adapted from RemoveSVNPaths, to remove arbitrary paths (BSH)
% 31.05.13 Allow empty string for pathList to get current path.
%          Remove filesep prepended to match string, so that what this does
%          matches the comment.

% If no pathList was passed to the function we'll just grab the one from
% Matlab.
if (nargin < 1 || isempty(pathList))
    % Grab the path list.
    pathList = path;
end

% if no matchString was provided, return the path unmodified
if nargin < 2
    newPathList = pathList;
    return;
end

try
    % We do the matching path removal in a try-catch block, because some of
    % the functions used inside this block are not available in Matlab-5
    % and GNU/Octave. Our catch - block provides fail-safe behaviour for
    % that case.
    
    % Break the path list into individual path elements.
    if IsOctave
        pathElements = strsplit(pathList,pathsep,true);
    else
        pathElements = textscan(pathList, '%s', 'delimiter', pathsep);
        pathElements = pathElements{1}.';
    end
    
    % Look at each element from the path.  If it doesn't contain a matching
    % folder then we add it to the end of our new path list.
    isNotMatching = cellfun(@isempty,strfind(pathElements,matchString));
    pathElements = pathElements(isNotMatching);
    
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
