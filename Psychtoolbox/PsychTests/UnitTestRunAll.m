% this script runs all m files in PsychTests/UnitTests, each of which is
% assumed to be a unittest that is to return true upon success or false
% otherwise

clear all; close all;

% get all scripts from subfolder (note that the FileFromFolder function
% can't be broken for this to work of course... :P)
[utests,ntests] = FileFromFolder(fullfile(PsychtoolboxRoot,'PsychTests/UnitTests'),'silent','m');
% delete contents file as its not executable
qContents = strcmp('Contents',{utests.fname});
utests(qContents) = [];
ntests = ntests-sum(qContents);

if ntests<1
    error('No Unit Tests found');
end

% build function handles for each
funs = {utests.fname};
funs = cellfun(@str2func,funs,'UniformOutput',false);

% run each
nfailed = 0;
fprintf('Running UnitTests\n');
for p=1:ntests
    fprintf('=========================================\n');
    fprintf('Running: PsychTests/UnitTests/%s:\n\n',utests(p).name);
    
    st = funs{p}();
    
    fprintf('\n');
    if st
        fprintf('PsychTests/UnitTests/%s: succes\n',utests(p).name);
    else
        fprintf('PsychTests/UnitTests/%s: failed\n',utests(p).name);
        nfailed = nfailed+1;
        failed(nfailed) = utests(p);
    end
    fprintf('=========================================\n');
    fprintf('\n');
end

fprintf('Summary:\n');
fprintf('%d tests succeeded, %d failed\n',ntests-nfailed,nfailed);
if nfailed>0
    fprintf('Tests failed:\n');
    fprintf('  %s\n',failed.fname);
end
