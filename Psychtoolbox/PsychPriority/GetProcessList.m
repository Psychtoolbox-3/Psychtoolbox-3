function processStructList=GetProcessList
% processStructList=GetProcessList
%
% OS X: ___________________________________________________________________
%
% Return an array of structures describing processes running on the system.
% GetProcessList ultimately relies on the Unix ps command.
%
% OS 9: ___________________________________________________________________
%
% GetProcessList does not exist in OS 9. 
%
% WINDOWS: ________________________________________________________________
% 
% GetProcessList does not exist in Windows.
% 
% _________________________________________________________________________
%
% see also: GetProcessDescriptorFromCommandName, GetRawProcessList

% HISTORY
%
% 12/09/03      awi     Wrote it.
% 07/15/04      awi     Added AssertOSX call, improved documention.    

% This script runs only on OS X.  It is used on OS X by Rush and Priority 
% which are platform neutral.  
AssertOSX;

%specify the names of fields which we will convert into integer doubles
%from the char reprenstations returned by ps.
numericFields={ 'uid', 'pid', 'ppid', 'cpu', 'pri', 'ni' ,'vsz','rss'};

% get an array of 
rawPlist=GetRawProcessList;
lineBreakCharacter=10;          %ascii newline
breakLocations=find(rawPlist==lineBreakCharacter);
lineStartIndices=[1 breakLocations(1:end-1)];
lineEndIndices=breakLocations;
lineArray={};
for i=1:length(lineStartIndices)
    lineArray{i}=rawPlist(lineStartIndices(i):lineEndIndices(i));
end
processLines=lineArray(2:end);
fieldNamesString=lineArray{1};

%get the field names from the ps column header
remainder=fieldNamesString;
fieldNames={};
gogo=1;
while gogo
    [tempFieldName, remainder]=strtok(remainder);
    if isempty(tempFieldName)
        gogo=0;
    else
        fieldNames={fieldNames{:} tempFieldName}; 
    end
end

%use lowercase field names so we dont' have to use the shift of caps lock
%keys when programming.
lcFieldNames=lower(fieldNames);

%fill fields of each index in the structure array values pulled from the
%corresponding line of the ps output.
for i=1:length(processLines)
    remainder=processLines{i};
    for j = 1:length(lcFieldNames)
        [token, remainder]=strtok(remainder);
        eval(['processStructList(i).' lcFieldNames{j} '=token;']);
    end
end

%convert the numeric fields into integer doubles from chars
for i = 1:length(numericFields)             %iterate over the field names
    if isfield(processStructList, numericFields{i})
        for j=1:length(processStructList)   %iterate over the array elements
            processStructList(j)=setfield(processStructList(j), numericFields{i}, sscanf(getfield(processStructList(j), numericFields{i}), '%d'));
        end %for
    end %if
 end %for




    
    





    
    
    
