function newRect = CenterRectOnPointd(rect,x,y)
% newRect = CenterRectOnPointd(rect,x,y)
% 
% CenterRectOnPointd offsets the rect to center it around an x and y position.
% The center of the rect is reported by RectCenterd. This is the same as
% CenterRectOnPoint, just that it returns floating point results instead of
% results that are rounded to integral pixel positions.
%
% rect can be a rect array, and/or x and/or y can be vectors instead of
% scalars.
%
% See also PsychRects/Contents, RectCenter, RectCenterd

% 9/13/99      Allen Ingling wrote it
% 10/6/99  dgp With permission from Allen, changed name from "PinRect", which 
%              conflicts with Apple's very different PinRect QuickDraw function, 
%              to CenterRectOnPoint. Simplified code.
% 5/18/08  mk  Vectorized, and uses RectCenterd for non-integral results.
% 7/27/15  dcn Fully support all vector/scalar combinations of rect x and y
%              inputs
% 7/29/15  dcn Now handles 4xM and Mx4 inputs

[cX,cY] = RectCenterd(rect);
newRect = OffsetRect(rect,x(:)-cX(:),y(:)-cY(:));
