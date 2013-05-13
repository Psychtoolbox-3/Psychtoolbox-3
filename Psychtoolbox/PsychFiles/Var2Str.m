function str = Var2Str(varargin)
% str = Var2Str(in,name)
%
% Takes variable IN and creates a string representation of it that would
% return the original variable when fed to eval(). NAME is the name of the variable
% that will be printed in this string.
% Can process any (combination of) MATLAB built-in datatype
%
% examples:
%   Var2Str({4,7,'test',@exp})
%   ans =
%       {4,7,'test',@exp}
%
%   var.field1.field2 = {logical(1),'yar',int32(43),[3 6 1 2],@isempty};
%   Var2Str(var,'var2')
%   ans =
%       var2.field1.field2 = {true,'yar',int32(43),[3,6,1,2],@isempty};

% DN 2008-01    Wrote it.
% DN 2008-07-30 Added support for function handles
% DN 2008-07-31 Added checking of MATLAB version where needed
% DN 2008-08-06 Lessened number of output lines needed to represent
%               variable in some cases
% DN 2008-08-12 Added wronginputhandler() for easy creation of specific
%               error messages on wrong input + added sparse support
% DN 2011-06-07 Simple variables while having only one input argument
%               didn't actually work....
% DN 2011-06-08 reworked 2D+ engine and added a dispatcher to remove code
%               duplication. Also put addition of LHS to expression in one
%               common place.
% DN 2011-06-09 Added handling of empty and/or fieldless structs

% TODO: make extensible by userdefined object parser for user specified
% datatypes - this will make this function complete

% make string of values in input, output is a cell per entry
strc = dispatcher(varargin{:});

if iscell(strc)
    % align equals-signs (autistic)
    indices     = strfind(strc,'=');
    indices     = cellfun(@(x)x(1),indices);
    MaxIndex    = max(indices);
    
    for p = 1:length(strc)
        i       = indices(p);
        dif     = MaxIndex - i;
        strc{p} = [strc{p}(1:i-1) repmat(' ',1,dif) strc{p}(i:end)];
    end
    
    % make string out of the stringcells
    str = [strc{:}];
else
    str = strc;
end


% dispatch string representation creation of datatype to handler for that
% datatype
function str = dispatcher(varargin)
if isnumeric(varargin{1}) || islogical(varargin{1})
    str = numeric2str(varargin{:});
elseif ischar(varargin{1})
    str = str2str(varargin{:});
elseif iscell(varargin{1})
    str = cell2str(varargin{:});
elseif isstruct(varargin{1})
    str = struct2str(varargin{:});
elseif isa(varargin{1},'function_handle')
    str = funchand2str(varargin{1});
else
    wronginputhandler(varargin{:});
end

% add LHS if needed
if ischar(str)
    if nargin==2
        str = {[varargin{2} ' = ' str ';' char(10)]};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions to deal with the different datatypes

% this function handles the unsupported data types
function wronginputhandler(input,name)
if nargin == 1
    at = '';
else
    at = [char(10) 'Error at ' name];
end
switch class(input)
    case 'inline'
        error('Inline functions not supported, they are deprecated.\nPlease see ''help function_handle'' for information on an alternative.%s',at)
    otherwise
        error('input of type %s is not supported.%s',class(input),at)
end

% this function handles numerical and logical data types
function str = numeric2str(input,name)

psychassert(isnumeric(input)||islogical(input),'numeric2str: Input must be numerical or logical, not %s',class(input))

if numel(input)>4 && isscalar(unique(input))
    % special case, all same value
    str = constant2str(input);
elseif isempty(input) || ndims(input)<=2
    if ismember(class(input),{'double','logical'})
        str = mat2str(input,17);
    else
        % for any non-double datatype, preserve type - double is default
        str = mat2str(input,17,'class');
    end
    str = regexprep(str,'\s+',' ');
else
    psychassert(nargin==2,'input argument name must be defined if processing 2D+ matrix');
    str = mat2strhd(input,name);
end

% this function handles the char data type -> strings
function str = str2str(input,name)

psychassert(ischar(input),'str2str: Input must be char, not %s',class(input))

if numel(input)>4 && isscalar(unique(input))
    % special case, all same value
    str = constant2str(input);
elseif isempty(input) || ndims(input)<=2
    str = mat2str(input,17);
else
    psychassert(nargin==2,'input argument name must be defined if processing 2D+ string matrix');
    str = mat2strhd(input,name);
end

% this one for cells
function str = cell2str(input,name)

psychassert(iscell(input),'cell2str: Input must be cell, not %s',class(input));

qstruct     = IsACell(input,@isstruct);         % recursively check if there is any struct in the cell
q2dplus     = IsACell(input,@(x)ndims(x)>2);    % recursively check if there is any element in the cell extending over more than 2 dimensions


if isempty(input)
    if ndims(input)==2 && all(size(input)==0)
        str = '{}';
    else
        s   = size(input);
        str = ['cell(' regexprep(mat2str(s,17),'\s+',',') ')'];
    end
elseif ~qstruct && ~q2dplus
    [nrow ncol] = size(input);
    str         = '{';
    % process cell per element
    for p = 1:nrow
        for q=1:ncol
            if isstruct(input{p,q})
                error('structs should not be processed here in any circumstance');
            else
                str = [str dispatcher(input{p,q})];
            end
            
            if q~=ncol
                str = [str ','];
            end
        end
        if p~=nrow
            str = [str ';'];
        end
    end
    str = [str '}'];
else
    psychassert(nargin==2,'input argument name must be defined if processing 2D+ cell or cell containing structs');
    str = cell2strhd(input,name);
end

% and this one for structs
function strc = struct2str(input,name)

if ~isstruct(input)
    error('Input is not struct')
end

%%%%%
fields  = fieldnames(input);
strc    = [];

if isempty(fields)
    % handle case of no fields
    strc = 'struct()';
    
    if ~isscalar(input)
        sizestr = regexprep(mat2str(size(input),17),'\s+',' ');
        strc    = ['repmat(' strc ',' sizestr ')'];
    end
elseif all(arrayfun(@(x) all(structfun(@(y) ndims(y)==2 && all(size(y)==0) && isa(y,'double'),x)),input))
    % handle case of all fields default-empty (0x0 double)
    nf = length(fields);
    fields = MergeCell('''',fields,''',[]');
    strc = cell2mat(['struct(' Interleave(fields,repmat(',',1,nf-1)) ')']);
    
    if ~isscalar(input)
        sizestr = regexprep(mat2str(size(input),17),'\s+',' ');
        strc    = ['repmat(' strc ',' sizestr ')'];
    end
else
    % struct has actual data, process
    psychassert(nargin==2,'input argument name must be defined if processing a struct');
    
    if isscalar(input)
        fields = fieldnames(input);
        for r=1:length(fields)
            wvar = input.(fields{r});
            namesuff = ['.' fields{r}];
            
            strc = [strc; dispatcher(wvar,[name namesuff])];
        end
    else
        strc = struct2strnonscalar(input,name);
    end
end



function str = funchand2str(input)

psychassert(isa(input,'function_handle'),'funchand2str: Input must be a function handle, not %s',class(input));

str = func2str(input);
if str(1)~='@'
    str = ['@' str];
end


%%%% function for special case of constant array
function str = constant2str(in)
item = unique(in);
s    = size(in);

sizestr = regexprep(mat2str(s,17),'\s+',' ');
itemstr = dispatcher(item);

if islogical(item)
    str = [itemstr '(' strrep(sizestr(2:end-1),' ',',') ')'];
else
    str = ['repmat(' itemstr ',' sizestr ')'];
end
return;



%%%% HD functions for variables of more than 2 non-singleton dimensions

function strc = mat2strhd(in,name)

s = size(in);
% unwrap all higher dimensions into a 2D mat, e.g., make a 4x2x3 into a
% 12x2 where every four rows contain one element from the third dimension
in=permute(in,[1,3:numel(s),2]);
in=reshape(in,[],s(2));
% wrap in cell per higher-dimension element
in=num2cellStrided(in,[s(1:2)]);

% prepare output indices
idxs    = size2idxs(s,2);
fmt     = [name '(:,:' repmat(',%d',[1,length(s)-2]) ')'];

% dispatch each 2D to correct interpreter
strc = cellfun(@(x,y)  dispatcher(x,[name sprintf(fmt,y)]),in,num2cell(idxs,2));

function strc = cell2strhd(in,name)

strc = [];
siz = size(in);

qDontProcess = cellfun(@(x) ndims(x)==2 && all(size(x)==0) && isa(x,'double'),in);
for p=1:numel(in)
    if qDontProcess(p)
        continue;
    end

    [idx{1:numel(siz)}] = ind2sub(siz,p);
    namesuff = ['{' Interleave([idx{:}],repmat(',',1,length(idx)-1)) '}'];
    
    strc = [strc; dispatcher(in{p},[name namesuff])];
end

function strc = struct2strnonscalar(in,name)

strc = [];
siz = size(in);

qDontProcess = arrayfun(@(x) all(structfun(@(y) ndims(y)==2 && all(size(y)==0) && isa(y,'double'),x)),in);
for p=1:numel(in)
    if qDontProcess(p)
        continue;
    end
    
    [idx{1:numel(siz)}] = ind2sub(siz,p);
    namesuff = ['(' Interleave([idx{:}],repmat(',',1,length(idx)-1)) ')'];
    
    strc = [strc; struct2str(in(p),[name namesuff])];
end


%%% other helpers
function res = size2idxs(siz,noff)
if nargin==1
    noff = 0;
end
narg=numel(siz);
n=narg-noff;

if n
    arg=cell(n,1);
    x=cell(n,1);
    for	i=1:n
        arg{i}=1:siz(i+noff);
    end
else
    res = [];
    return;
end

if n > 1
    [x{1:n,1}]=ndgrid(arg{1:end});
    res=reshape(cat(n+1,x{:}),[],n);
else
    res=arg{:}.';
end

function res = num2cellStrided(in,stride)
% the function is only tested 2D for now, don't want to think about higher
% dims, but it might just work equally fine
siz     = size(in);
ndim    = numel(siz);
ncell   = siz./stride;
assert(~any(mod(ncell,1)),'size needs to be a multiple of stride for each dimension');

res = reshape( in , Interleave(stride,ncell) );
res = permute( res, InterLeave(1:ndim,ndim+1:ndim*2) );

res = squeeze(num2cell(res,[1:ndim]));
