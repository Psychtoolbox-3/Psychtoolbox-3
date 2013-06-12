function directory=CalDataFolder(forceDemo,calFileName,calDir)
% directory=CalDataFolder([forceDemo],[calFileName],[calDir])
%
% Get the path to the CalData folder.
%
% If "forceDemo" is true (false by default), then force use of
% PsychCalDemoData.  Otherwise return location of PsychCalLocalData
% if it exists, and PsychCalDemoData if it PsychCalLocalData does
% not exist.
%
% If calFileName is passed, it checks if the specified file exists
% in the top level of the directory it located.  If so, it returns
% that directory.  If not, it searches the first level subdirectories
% of the top level directories and if the specfied file is in one of
% them, returns that.  Skips subdirs called 'xOld', 'Plots', and those
% that begin with '.'.
%
% calFileName can either have .mat postpended or not.
%
% If calDir is passed as a string, that is used as the directory, with
% the subdir searching as above.
%
% See also LoadCalFile, SaveCalFile.
%
% Denis Pelli 7/25/96
% Denis Pelli 2/28/98 change "CalDat" to "PsychCalData"
% 8/14/00  dhb  Add alternate name, change names. 
% 4/1/07   dhb  Fix subtle bug in error message when there are duplicate cal
%               folders on path. 
% 3/7/08   mpr  changed documentation to make it consistent (apparently
%               "forceDemo" used to be "alt"
% 4/2/13   dhb  Add calFileName and associated behavior.
% 6/2/13   dhb  Make this properly return subfolder containing calibration file
%               if .mat is not postpended.
% 6/10/13  dhb  Fix buglet introduced 6/2/13 -- need to handle empty calFileName (thanks to MS for
%               identifying the problem and the fix.

% Set forceDemo flag
if (nargin < 1 || isempty(forceDemo))
	forceDemo = 0;
end
if (nargin < 2 || isempty(calFileName))
    calFileName = [];
end

% Postpend .mat if necessary
if (~isempty(calFileName) && ~strcmp(calFileName(end-3:end),'.mat'))
    calFileName = [calFileName '.mat'];
end

% If dir is passed we just use that.  Otherwise
% do our thing.
if (nargin < 3 || isempty(calDir)) 
    name='PsychCalLocalData';
    alternateName ='PsychCalDemoData';
    
    % Find name.  If not there, find alternate
    if (~forceDemo)
        directory = FindFolder(name);
        duplicateMsgName = name;
    else
        directory = [];
    end
    if isempty(directory)
        directory=FindFolder(alternateName);
        duplicateMsgName = alternateName;
    end
    
    % If both finds fail, print out error message.  This
    % should never happen as we put 'PsychCalDemoData' in
    % the toolbox distribution.
    if isempty(directory)
        error(['Can''t find any ''' name ''' or ''' alternateName '''folders in the Matlab path.']);
    end
    
    % If we found multiple copies of a calibration folder, we complain.
    % This also should never happen.
    if size(directory,1)>1
        for i=1:size(directory,1)
            disp(['DUPLICATE: ''' deblank(directory(i,:)) '''']);
        end
        error(['Found more than one ''' duplicateMsgName ''' folder in the Matlab path.']);
    end
else
    directory = calDir;
end

% If we were passed a calibration file name, figure out whether it exists in
% the found directory, or a subfolder of the directory.
if (~isempty(calFileName))
    curDir = pwd;
    cd(directory);
    
    % Is the file we're after here?  If so, we're good
    % and we just return after restoring directory.
    % Otherwise, we look in any subdirectories other
    % than those beginning with 'xOld' or 'Plot' and
    % see if our baby is there.
    foundOne = 0;
    dirRet = dir([calFileName]);
    if (~isempty(dirRet))
        foundOne = 1;
    end
    allDirs = dir;
    for i = 1:length(allDirs)
        if (allDirs(i).isdir && ~strcmp(allDirs(i).name,'xOld') && ~strcmp(allDirs(i).name,'Plots') && allDirs(i).name(1) ~= '.')
            cd(allDirs(i).name);
            dirRet = dir(calFileName);
            cd('..');
            if (~isempty(dirRet))
                if (foundOne)
                    error('More than one copy of the desired calibration file found');
                end
                directory = [directory allDirs(i).name filesep];
                foundOne = 1;
            end
        end
    end
    
    cd(curDir);
end
    
