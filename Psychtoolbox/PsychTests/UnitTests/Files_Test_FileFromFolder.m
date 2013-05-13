function success = Files_Test_FileFromFolder

success = true;

% create some files in a directory, for testing
newdir = fullfile(PsychtoolboxRoot,'PsychTests','UnitTests','FileTests');
mkdir(newdir);
fid=fopen(fullfile(newdir,'test1.txt'),'a'); fclose(fid);
fid=fopen(fullfile(newdir,'test2.txt'),'a'); fclose(fid);
fid=fopen(fullfile(newdir,'test1.gif'),'a'); fclose(fid);
fid=fopen(fullfile(newdir,'test3_noext'),'a'); fclose(fid);
mkdir(fullfile(newdir,'testdir'));

try
    % test: FileFromFolder
    [f,nf]=FileFromFolder(newdir);
    if nf~=4
        success = false;
        fprintf('FileFromFolder didn''t get right number of files, got:\n');
        fprintf('  %s\n',f.fname);
    elseif any([f.isdir])
        success = false;
        fprintf('FileFromFolder got directories\n');
    end
    clear f; clear nf
    
    [f,nf]=FileFromFolder(newdir,'','gif');
    if nf~=1 || ~(strcmp(f.name,'test1.gif') && strcmp(f.fname,'test1') && strcmp(f.ext,'gif'))
        success = false;
        fprintf('FileFromFolder didn''t get the right file when using file extension mask\n');
    end
    clear f; clear nf
    
    [f,nf]=FileFromFolder(newdir,'','');
    if nf~=1 || ~(strcmp(f.fname,'test3_noext') && isempty(f.ext))
        success = false;
        fprintf('FileFromFolder didn''t get the right file when extension mask to find files without extension\n');
    end
    clear f; clear nf
    
    % test: FolderFromFolder
    [f,nf]=FolderFromFolder(newdir);
    if nf~=1 || ~f.isdir || ~strcmp(f.name,'testdir')
        success = false;
        fprintf('FolderFromFolder didn''t get the right folder\n');
    end
    
catch
    success = false;
    fprintf('error ocurred: "%s"\n',lasterr);
    
end

% cleanup
if IsOctave, OLD_VAL = confirm_recursive_rmdir (false); end
rmdir(newdir,'s');
if IsOctave, confirm_recursive_rmdir(OLD_VAL); end
