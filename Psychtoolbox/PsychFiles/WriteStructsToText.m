function WriteStructsToText(filename,theStructs)
% WriteStructsToText(filename,theStructs)
%
% Write a tab delimited text file.  The first row should
% contain the field names for a structure.  Each following
% row contains the data for one instance of that structure.
%
% This routine writes each element of the structure array as a row.
% Only numeric and string field values are supported.
%
% If the filename is a string, this routine handles opening and
% closing of the file.  Otherwise, this routine assumes
% that what was passed is a valid fid and uses it as such,
% leaving the calling routine to do the opening and closing.
%
% 06/16/03 dhb  Wrote it.
% 07/01/03 dhb  Support string fields too.
% 08/12/06 dhb  Call filetype only if OS9.
% 07/07/13 dhb  Add handling of non-string filename as fid.

% Open the file
if (isstr(filename))
    fid = fopen(filename,'wt');
else
    fid = filename;
end
if (fid == -1)
    error('Error opening file or invalid fid passed');
end

% Get the fieldnames
theFields = fieldnames(theStructs(1));
nFields = length(theFields);
for i = 1:nFields
	fprintf(fid,'%s',theFields{i});
	if (i < nFields)
		fprintf(fid,'\t');
	else
		fprintf(fid,'\n');
	end
end

% Now write each struct's data as a line
nStructs = length(theStructs);
for j = 1:nStructs
	for i = 1:nFields	
		if (ischar(getfield(theStructs(j),theFields{i})))
			fprintf(fid,'%s',getfield(theStructs(j),theFields{i}));
		else
			fprintf(fid,'%g',getfield(theStructs(j),theFields{i}));
		end
		if (i < nFields)
			fprintf(fid,'\t');
		else
			fprintf(fid,'\n');
		end
	end
end

% Close the file.
if (isstr(filename))
    fclose(fid);
end
