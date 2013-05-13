function isEmpty = IsEmptyRect(a)
% isEmpty = IsEmptyRect(a)
% 
% Returns 1 if the rect "a" is empty, otherwise returns 0.
% "A rectangle is considered empty if the bottom coordinate is less than 
% or equal to the top coordinate or if the right coordinate is less than 
% or equal to the left." (Apple, Inside Mac)
% Also see PsychRects.

% 7/10/96 dgp  Wrote it.
% 8/5/96 dgp check rect size.

if nargin~=1
	error('Usage:  isEmpty = IsEmptyRect(a)');
end
if size(a,2)~=4
	error('Wrong size rect argument. Usage:  isEmpty = IsEmptyRect(a)');
end
isEmpty=RectWidth(a)<=0 || RectHeight(a)<=0;
