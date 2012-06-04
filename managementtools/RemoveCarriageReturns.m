function RemoveCarriageReturns(folder,qRecursive,qDryRun, exts_filter)
% finds all 13s and 10s (decimal char number)
% throws out all 13s, but if not followed by 10 it replaces it by 10
%
% make sure to do a dryrun first and set the fourth argument such that no
% files you dont want to touch get touched. e.g.: {'m','c'} to process only
% *.m and *.c files in the tree


if nargin<4
    exts_filter = {};
else
    exts_filter = {exts_filter};
end

if nargin<3
    qDryRun = true;
end

if nargin<2
    qRecursive = false;
end

if qRecursive
    [folds,nfolds] = FolderFromFolder(folder,'ssilent');
    for p=1:nfolds
        RemoveCarriageReturns(fullfile(folder,folds(p).name),qRecursive,qDryRun, exts_filter{:});
    end
end

[files,nfiles] = FileFromFolder(folder,'ssilent',exts_filter{:});

if nfiles==0
    return;
end
   
for p=1:nfiles
    % open it
    fname = fullfile(folder,files(p).name);
    fid = fopen(fname,'rt');
    dat = fread(fid);
    fclose(fid);
    
    % do replace
    newdat = DoReplace(dat);
    
    % save if changed
    if ~isequal(dat,newdat)
        fprintf('%s\n',files(p).name);
        
        if ~qDryRun
            fid = fopen([folder filesep files(p).name],'wt');
            fwrite(fid,newdat);
            fclose(fid);
        end
    end
end



function dat = DoReplace(dat)

% find all 13s
q13 = dat==13;
% find all 10s
q10 = dat==10;

% find all 13s followed by 10s
q = q13(1:end-1) & q10(2:end);

% remove them
dat(q) = [];

% all other 13s are orphans, replace them by 10s
dat(dat==13) = 10;