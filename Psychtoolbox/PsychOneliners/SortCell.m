function [y,ndx] = SortCell(x, dim)
% SortCell    Sort a cell array in ascending order.
%
% Description: SortCell sorts the input cell array according to the
%   dimensions (columns) specified by the user.
%
% Usage: [Y,NDX] = SortCell(X, DIM)
%
% Input:
%	 X: the cell array to be sorted.
%  DIM: (optional) one or more column numbers. Simply an array of one or
%       more column numbers.  The first number is the primary column on
%       which to sort. Extra column numbers may be supplied if secondary
%       sorting is required. The default value is 1, if no dimension
%       array is supplied.
%
% Output:
%     Y: the sorted cell array.
%   NDX: indices such that Y = X(NDX,:).
%
% Example:    Y = SortCell(X, [3 2])
%
% Note that this function has only been tested on mixed cell arrays
% containing character strings and numeric values.

%   Copyright 2007 Jeff Jackson (Ocean Sciences, DFO Canada)
%   Creation Date: Jan. 24, 2007
%   Last Updated:  Jan. 25, 2007
%   2008 DN v1.1:  Added support for numerical datatypes other than char
%                  and double
%   2016 DN v1.2:  Now also outputs index into input to create output.


% Check input arguments
if nargin == 0
    error('no input arguments were supplied.  at least one is expected.');
elseif nargin == 1
    dim = 1;
elseif nargin == 2
    if ~isnumeric(dim)
        error('the second input argument is not numeric.  at least one numeric value is expected.');
    end
else
    error('too many input arguments supplied.  only two are allowed.');
end
if ~iscell(x)
    error('the first input argument is not a cell array.  a cell array is expected.');
end

% prepare output
ndx = [1:size(x,1)].';

% Now find out if the cell array is being sorted on more than one column.
% If it is then use recursion to call the SortCell function again to sort
% the less important columns first. Repeat calls to SortCell until only one
% column is left to be sorted. Then return the sorted cell array to the
% calling function to continue with the higher priority sorting.
ndim = length(dim);
if ndim > 1
    col     = dim(2:end);
    [x,xi]  = SortCell(x, col);
    ndx     = ndx(xi);
end

% Get the dimensions of the input cell array.
nrows   = size(x,1);

% Retrieve the primary dimension (column) to be sorted.
col     = dim(1);

% Place the cells for this column in variable 'b'.
b       = x(:,col);

% Check each cell in cell array 'b' to see if it contains either a
% character string or numeric value.
qchar   = cellfun(@(x)isa(x,'char') , b);
classes = cellfun(@class            , b,'UniformOutput',false);

% Check to see if cell array 'b' contained only character strings.
% If cell array b contains character strings then do nothing because
% no further content handling is required.
if sum(qchar) == nrows
    % Check to see if cell array 'b' contained only numeric values of
    % the same type.
elseif length(unique(classes))==1 && ismember(unique(classes),{'logical','single','double','float','int8','uint8','int16','uint16','int32','uint32','int64','uint64'})
    % If the cells in cell array 'b' contain numeric values retrieve the cell
    % Contents and change 'b' to a numeric array.
    b = [b{:}];
else
    error('This column (%d) is mixed so sorting cannot be completed.',dim(1));
end

% Sort the current array and return the new index.
[~,xi] = sort(b);

% Using the index from the sorted array, update the input cell array and
% return it.
y   = x(xi,:);
ndx = ndx(xi,:);
