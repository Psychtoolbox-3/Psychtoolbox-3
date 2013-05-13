function [A, tf] = Replace(A, S1, S2)
% Replace - Replace Elements
%   B = Replace(A,S1,S2) returns a matrix B in which the elements in A that 
%   are in S1 are Replaced by those in S2. In general, S1 and S2 should have
%   an equal number of elements. If S2 has one element, it is expanded to
%   match the size of S1. Examples:
%      Replace([1 1 2 3 4 4],[1 3],[0 99]) % ->  [ 0 0 2 99 4 4]
%      Replace(1:10,[3 5 6 8],NaN) % ->  [ 1 2 NaN 4 NaN NaN 7 NaN 9 10]
%      Replace([1 NaN Inf 8 99],[NaN Inf 99],[12 13 14]) % -> [1 12 13 8 14]
%
%   [B, TF] = Replace(A,S1,S2) also returns a logical vector TF of the same
%   size as A. TF is true for those elements that are replaced.
%
%   A and S1 can be cell arrays of strings. In that case S2 should be a
%   cell array as well but can contain mixed types. Example:
%      Replace({'aa' 'b' 'c' 'a'},{'a' 'b'}, {'xx' 2}) %-> {'aa' [2] 'c' 'xx'}
%
%   If S2 is empty, the elements of A that are in S1
%   are removed. Examples:
%      Replace(1:5,[2 4],[]) % -> [1 3 5]
%      Replace({'aa' 'a' 'b' 'aa' 'c'},{'aa','c'},{}) % -> {'a', 'b'}
%
%   See also FIND, STRREP, REGEXPREP, ISMEMBER

% for Matlab R13
% version 1.4 (dec 2006)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.0 (feb 2006) created
% 1.1 (feb 2006) fixed bug when NaNs were to be removed 
% 1.2 (feb 2006) fixed again bug with NaNs
% 1.3 (oct 2006) fixed error when using matrices
% 1.4 (dec 2006) added additional outputs of TF and LOC
% 1.5 (apr 2008) DN: optimized for R2007b

error(nargchk(3,3,nargin)) ;

% all three inputs should be cell arrays or numerical arrays
if ~isequal(iscell(A), iscell(S1), iscell(S2)),
    error('The arguments should be all cell arrays or not.') ;
end

if iscell(A),
    % if they are cell, they should be character arrays
    if ~iscellstr(A),
        error('A should be a cell array of strings.') ;
    end
    if ~iscellstr(S1),
        error('S1 should be a cell array of strings.') ;
    end
end

if ~isempty(S2),
    if numel(S2)==1,
        % single element expansion
        S2 = repmat(S2,size(S1)) ;
    elseif numel(S1) ~= numel(S2),
        error('The number of elements in S1 and S2 do not match ') ;
    end
end

% the engine
[tf,loc] = ismember(A(:),S1(:)) ;

if nargout>1,
    tf = reshape(tf,size(A)) ;
end

if any(tf),
    if isempty(S2),
        A(tf) = [] ;
    else
        A(tf) = S2(loc(tf)) ;
    end
end

% special treatment for nans if necessary
if ~iscell(S1),
    % only for non-cell arrays
    qsn = isnan(S1(:)) ;
    if any(qsn),
        qa = isnan(A(:)) ;
        if any(qa),            
            if isempty(S2),
                A(qa) = [] ;
            else
                i = find(qsn,1,'first') ;            
                A(qa) = S2(i) ;
            end
        end
    end
end
