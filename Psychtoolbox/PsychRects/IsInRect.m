function inside = IsInRect(x,y,rect)
% inside = IsInRect(x,y,rect)
%
% Is location x,y inside the rect?
%
% Also see PsychRects.

% 3/5/97  dhb  Wrote it.

if (x >= rect(RectLeft) && x <= rect(RectRight) && ...
		y >= rect(RectTop) && y <= rect(RectBottom) )
	inside = 1;
else
	inside = 0;
end
