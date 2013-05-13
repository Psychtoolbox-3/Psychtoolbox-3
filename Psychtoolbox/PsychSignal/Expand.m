function B=Expand(A,horizontalFactor,verticalFactor)
% B=Expand(A,horizontalFactor,[verticalFactor])
%
% Expands the ND matrix A by cell replication, and returns the result.
% If the vertical scale factor is omitted, it is assumed to be 
% the same as the horizontal. Note that the horizontal-before-vertical
% ordering of arguments is consistent with image processing, but contrary 
% to Matlab's rows-before-columns convention.
%
% We use "Tony's Trick" to replicate a vector, as explained
% in MathWorks Matlab Technote 1109, section 4.
%
% Also see ScaleRect.m

% Denis Pelli 5/27/96, 6/14/96, 7/6/96
% 7/24/02 dgp Support an arbitrary number of dimensions.
% 13/06/12 DN Redid internals for significant speedup.

if nargin<2 || nargin>3
	error('Usage: A=Expand(A,horizontalFactor,[verticalFactor])');
end
if nargin==2
	verticalFactor=horizontalFactor;
end

psychassert(round(verticalFactor)  ==verticalFactor   && verticalFactor>=1 && ... 
            round(horizontalFactor)==horizontalFactor && horizontalFactor>=1, ...
        	'Expand only supports positive integer factors.');
psychassert(~isempty(A),'Can''t expand an empty matrix');


% Generate row copying instructions index.
inds                            = 1:size(A,1);
rowCopyingInstructionsIndex     = inds(ones(verticalFactor,1),:);
rowCopyingInstructionsIndex     = rowCopyingInstructionsIndex(:);

% Generate column copying instructions index.
inds                            = 1:size(A,2);
columnCopyingInstructionsIndex  = inds(ones(horizontalFactor,1),:);
columnCopyingInstructionsIndex  = columnCopyingInstructionsIndex(:)';

% The following code uses Matlab's matrix indexing quirks to magnify the
% matrix.  It is easier to understand how it works by looking at a specific
% example:
% 
% >> n = [1 2; 3 4] % Matlab, please give me a matrix with four elements.
%
% n =
% 
%      1     2
%      3     4
% 
% >> % Matlab, please generate a new matrix by using the provided copying
% >> % instructions index.  My copying instructions index says that you
% >> % should print the first column twice, then print the second column
% >> % twice.  Thanks.
% >> m = n(:, [1 1 2 2])
% 
% m =
% 
%      1     1     2     2
%      3     3     4     4
%
% >> % Matlab, please generate a new matrix by using the provided copying
% >> % instructions index.  My copying instructions index says that you
% >> % should print the first row twice, then print the second row
% >> % twice.  Thanks.
% >> m = n([1 1 2 2], :)
% 
% m =
% 
%      1     2
%      1     2
%      3     4
%      3     4
%

if ndims(A)>3
    colon = {':'};
    B = A(rowCopyingInstructionsIndex, columnCopyingInstructionsIndex, colon{ones(1,ndims(A)-2)});
else
    B = A(rowCopyingInstructionsIndex, columnCopyingInstructionsIndex,:);
end
