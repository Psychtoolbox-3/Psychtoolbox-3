function [y,ia,ic] = UniqueCell(x)
% UniqueCell    Get unique rows from cell array
%
% Description: like unique(x,'rows'), but supports cell arrays, and these
%   cell arrays can have columns with different types
%
% Usage: [Y,IA,IC] = UniqueCell(X)
%
% Input:
%	 X: the cell array to be processed.
%
% Output:
%     Y: contains the same values as in A but with no repetitions. Y will
%        be sorted.
%    IA: index vector IA such that C = A(IA,:).
%    IC: index vector IC such that A = C(IC,:).
%
% Note that this function has only been tested on mixed cell arrays
% containing character strings and numeric values.

%   Copyright 2007 Diederick C. Niehorster (Lund University, Sweden)


% Check input arguments
if ~iscell(x)
    error('the first input argument is not a cell array.  a cell array is expected.');
end

% sort input, get indices
[sortX,indSortX] = SortCell(x,1:size(x,2));

% per column, find out where the values change
groupsSortX = false(size(sortX,1)-1,size(sortX,2));
for p=1:size(sortX,2)
    if iscellstr(sortX(:,p))
        groupsSortX(:,p) = ~strcmp(sortX(1:end-1,p),sortX(2:end,p));
    else
        groupsSortX(:,p) = cat(1,sortX{1:end-1,p}) ~= cat(1,sortX{2:end,p});
    end
end
groupsSortX = any(groupsSortX,2);
groupsSortX = [true; groupsSortX];  % First row is always a member of unique list.

% Extract Unique elements.
y = sortX(groupsSortX,:);           % Create unique list by indexing into sorted list.

% Find ia.
if nargout > 1
    ia = indSortX(groupsSortX);     % Find the indices of the sorted logical.
end

% Find ic.
if nargout == 3
    ic = cumsum(groupsSortX);       % Lists position, starting at 1.
    ic(indSortX) = ic;              % Re-reference indC to indexing of sortA.
end
