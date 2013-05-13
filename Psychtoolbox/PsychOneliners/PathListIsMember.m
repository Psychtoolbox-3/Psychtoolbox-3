function pathList = PathListIsMember(array,set)

% If no pathList set was passed to the function we'll just grab the one
% from Matlab.
if nargin < 2
    % Grab the path list.
    set = path;
end

% keep old copy of array in case the below fails (see catch)
arrayBackup = array;

try
    % We do the .svn path removal in a try-catch block, because some of the
    % functions used inside this block are not available in Matlab-5 and
    % GNU/Octave. Our catch - block provides fail-safe behaviour for that
    % case.
    
    % Break the path lists into individual path elements.
    if IsOctave
        array = strsplit(array,pathsep,true);
        set   = strsplit(set  ,pathsep,true);
    else
        array = textscan(array, '%s', 'delimiter', pathsep);
        array = array{1}.';
        set   = textscan(set  , '%s', 'delimiter', pathsep);
        set   = set{1}.';
    end
    
    % Look at each element from array and see if it is in set.  If not,
    % remove it from the pathlist
    qInSet = ismember(array,set);
    array = array(qInSet);
    
    if ~isempty(array)
        % generate new pathList
        array = [array; repmat({pathsep},1,length(array))];
        pathList  = [array{:}];
        
        % drop last separator
        pathList(end) = [];
    end
catch
    % Fallback behaviour: We fail-safe by simply returning the unmodified
    % array input. Nothing removed, but the whole beast is still
    % functional.
    pathList = arrayBackup;
end
