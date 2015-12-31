function isEmpty = IsEmptyRect(a)
% isEmpty = IsEmptyRect(a)
% 
% Returns true if the rect "a" is empty, otherwise returns false. If rect
% is a rect-array, returns as above for each rect in the array.
% "A rectangle is considered empty if the bottom coordinate is less than 
% or equal to the top coordinate or if the right coordinate is less than 
% or equal to the left." (Apple, Inside Mac)
%
% See also PsychRects/Contents.

% 7/10/96 dgp  Wrote it.
% 8/5/96  dgp  check rect size.
% 7/26/15 dcn  vectorized
% 7/29/15 dcn  Removed check for rect-array orientation, both are handled
%              fine.

if nargin~=1
	error('Usage:  isEmpty = IsEmptyRect(a)');
end
isEmpty=RectWidth(a)<=0 | RectHeight(a)<=0;
