function [width,height] = RectSize(rect)
% [width,height] = RectSize(rect)
%
% Returns the rect's width and height, or the width and height of each rect
% in the rect array.
%
% See also PsychRects/Contents, SizeOfRect.

% 10/10/2000 fwc wrote it based on PsychToolbox RectHeight/RectWidth
% 07/29/2015 dcn Removed check for rect-array orientation, both are handled
%                fine.

if nargin~=1
	error('Usage:  [width,height] = RectSize(rect)');
end

width = RectWidth(rect);
height = RectHeight(rect);
