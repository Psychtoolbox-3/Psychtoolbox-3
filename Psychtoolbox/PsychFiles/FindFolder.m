function directory=FindFolder(name)
% directory=FindFolder(name)
% Searches the Matlab 'path' for the path to the named folder.
% There may be no matches, one match, or multiple matches.
% Each unique match appears as a row in "directory".
% If there are no matches then "directory" will be an empty matrix.
% Matching ignores case.
% You should DEBLANK a row of "directory" before using it.
% 
% Also see Matlab's MKDIR, TEMPDIR, ISDIR, LOOKFOR, WHAT.
% Try HELP PsychFiles.

% 7/26/96  dgp Wrote it.
% 12/10/01 awi Set ignoreCase to 1 always.
% 4/13/02  dgp Updated to use Matlab's predefined separator symbols.

% Matlab predefines these:
% PATHSEP = path separator character
% FILESEP = directory separator character

ignoreCase=1;
paths=[path pathsep];
if ignoreCase
	n=lower(name);
	p=lower(paths);
else
	n=name;
	p=paths;
end
pathIndex=[1 1+findstr(pathsep,p)];
nameIndex=[findstr([filesep n filesep],p) findstr([filesep n pathsep],p)];
clear n p
if isempty(nameIndex)
	directory=[];
else
	nIndex=nameIndex(1);
	pIndex=pathIndex(max(find(pathIndex<=nIndex)));
	directory=[paths(pIndex:nIndex+length(name)) filesep];
	for i=2:length(nameIndex)
		nIndex=nameIndex(i);
		pIndex=pathIndex(max(find(pathIndex<=nIndex)));
		new=[paths(pIndex:nIndex+length(name)) filesep];
		unique=1;
		for j=1:size(directory,1)
			if streq(deblank(directory(j,:)),new)
				unique=0;
				break
			end
		end
		if unique
			directory=char(directory,new);
		end
	end
end

% FindFolder always succeeds. "directory" will have zero to many rows,
% each row containing a unique match. 

% What follows is optional code that you may want to add to your program, after
% it calls FindFolder, to give an error unless there was exactly one match.
if 0
	if isempty(directory)
		error(['Can''t find any ''' name ''' folder in the Matlab path.']);
	end
	if size(directory,1)>1
		for i=1:size(directory,1)
			disp(['DUPLICATE: ''' deblank(directory(i,:)) '''']);
		end
		error(['Found more than one ''' name ''' folder in the Matlab path.']);
	end
end
