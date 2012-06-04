function [file,nfile] = FileFromFolder(folder,mode,f_ext)
% [file,nfile] = FileFromFolder(folder,mode,ext)
%
% Returns struct with all files in directory FOLDER.
% MODE specifies whether an error is displayed when no directories are
% found (default). If MODE is 'silent', only a message will will be
% displayed in the command window, if 'ssilent', no notification will be
% produced at all. If left empty, default is implied.
% Ext is an optional filter on file extension. If specified, only files
% with the specified extension will be found. It can be a cell vector of
% strings for filtering on multiple extensions

% 2007 IH        Wrote it.
% 2007 IH&DN     Various additions
% 2008-08-06 DN  All file properties now in output struct
% 2009-02-14 DN  Now returns all files except '..' and '.', code
%                optimized
% 2010-05-26 DN  Got rid of for-loop, added optional filter on extension
% 2010-05-30 DN  Woops, some of the new changes break the function when no
%                files are found
% 2010-07-02 DN  Now supports filtering on multiple extensions
% 2010-07-12 DN  Fixed . at end of fname
% 2011-06-07 DN  Can now also filter for files with no extension
% 2012-06-04 DN  Now also have ssilent mode for no output at all

if nargin >= 2 && strcmp(mode,'silent')
    silent = 1;
elseif nargin >= 2 && strcmp(mode,'ssilent')
    silent = 2;
else
    silent = 0;
end


file        = dir(folder);
file        = file(~[file.isdir]);  % get rid of folders. This also skips '..' and '.', which are marked as dirs

if ~isempty(file)
    % get file name and extension
    [name,ext]  = cellfun(@SplitFName,{file.name},'UniformOutput',false);
    [file.fname]= name{:};
    [file.ext]  = ext{:};

    % if filter, use it
    if nargin >= 3
        q_ext   = ismember(ext,f_ext);
        file    = file(q_ext);
    end
end

nfile       = length(file);

if nfile==0
    if silent==1
        fprintf('FileFromFolder: No files found in: %s\n',folder);
        file = [];
    elseif ~silent
        error('FileFromFolder: No files found in: %s',folder);
    end
end



% helpers
function [name,ext] = SplitFName(name)
% Look for EXTENSION part
ind = find(name == '.', 1, 'last');

if isempty(ind)
    ext = '';
else
    ext = name(ind+1:end);
    name(ind:end) = [];
end
