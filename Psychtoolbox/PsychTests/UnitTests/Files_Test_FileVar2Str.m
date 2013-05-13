function success = Files_Test_FileVar2Str

if IsOctave
    fprintf('Var2Str is not supported on Octave as its mat2str doesn''t handle logical or sparse input and we depend on that (I don''t want to replicate functionality from matlab''s mat2str)\nVar2Str will work ok on Octave if you input numeric non-empty data, possibly contained in cells or structs.\n')
    success = false;
    return;
end

success = true;

try
    % create some simple and complicated variables, Var2Str should be able
    % to handle all of them.
    a = [3 4 5];
    b = [true false];
    c = char([116 101 115 116 32 13]);
    d = @sin;
    e = @(x) sin(x);    % anonymous function handle
    f = {1,2, [3 4 5], true, uint32([335321 32]), zeros([0, 3],'int16')};
    f{50,1,3,45} = {6,7,8,'a'};
    g.f1 = true;
    g(5).get = 'get';
    g(30).cell{1,4,62,2,5}.test(45,2,10).cell2{3,3,54,2,1} = {4,7,'test',@exp};
    % compression test (above we already that default-empty values are not
    % output, now test if fields of the same number (constant arrays) are
    % compressed
    h{1}  = {[6 6 6 6 6]};
    h{5}  = {int32([6 6 6 6 6])};
    h{9}  = false(5,7,3,2);
    h{13} = repmat('a',[3,4,23,2,7]);
    h{15} = zeros(30,30,5);
    h{18} = zeros(30,30,5,'uint16');
    h{19} = sprintf('\n\r');
    h{21} = speye(7);
    h{23} = [nan+1j*inf 0 0 0 0 inf-1j*inf 0 0 0 2.3-1j*pi];
    h{25} = sparse(h{23});
    h{27} = cell(0,3);
    h{29} = repmat(struct('f1',[],'f2',[]),[0 3]);
    h{30} = repmat(struct('f1',[],'f2',[]),[1 3 5 20 1]);
    h{31} = struct();
    h{32} = repmat(struct(),[3 4]);
    h{33} = repmat(struct(),[0 4]);
    
    % when evaling the output of Var2Str, we should get exactly the same
    % data back, lets do it
    success = Var2StrTester(a,'a') && success;
    success = Var2StrTester(b,'b') && success;
    success = Var2StrTester(c,'c') && success;
    success = Var2StrTester(d,'d') && success;
    success = Var2StrTester(f,'f') && success;
    success = Var2StrTester(g,'g') && success;
    success = Var2StrTester(h,'h') && success;
    
    % anonymous function handles won't compare equal, so need to deal with
    % this by hand
    str = Var2Str(e,'e2');
    eval(str);
    f1 = functions(e);
    f2 = functions(e2);
    if ~strcmp(f1.function,f2.function(end+[-length(f1.function)+1:0]))
        success = false;
        fprintf('Var2Str failed on variable e\n');
    end
    
    % check the simple one input version works as well
    if ~isequal(a,eval(Var2Str(a)))
        success = false;
        fprintf('Var2Str failed with 1 input version\n');
    end
    
    % and check LogVar (only check a trivial one and a complicated one as
    % LogVar is only a thin wrapper around Var2Str)
    % This fails when no write access to PTB tree
    success = LogVarTester(a,'a') && success;
    success = LogVarTester(h,'h') && success;
    
    
catch
    success = false;
    err = lasterror;
    fprintf('error ocurred: "%s"\nstack:\n',err.message);
    for p=1:length(err.stack)
        fprintf('Error in ==> %s at %d\n',err.stack(p).name,err.stack(p).line);
    end
end





function [success,var2,str] = Var2StrTester(var,name) %#ok<STOUT>
str = Var2Str(var,'var2');
eval(str);
if ~isequalwithequalnans(var,var2)
    success = false;
    fprintf('Var2Str failed on variable %s\n',name);
else
    success = true;
end

function [success,var2,str] = LogVarTester(var,name) %#ok<STOUT>
fname = LogVar(var,'var2',cd);
fid = fopen(fullfile(cd,fname),'rt');
str = fread(fid,inf,'*char');
fclose(fid);
delete(fullfile(cd,fname));
eval(str);
if ~isequalwithequalnans(var,var2)
    success = false;
    fprintf('LogVar failed on variable %s\n',name);
else
    success = true;
end
