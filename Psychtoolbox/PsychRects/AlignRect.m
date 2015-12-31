function rect=AlignRect(rect,fixedRect,side1,side2)
% rect=AlignRect(rect,fixedRect,side1,[side2])
%
% Moves rect to align its center/top/bottom/left/right with the
% corresponding edge(s)/point of fixedRect. Does "side1" and then "side2".
% The legal values for side1 and side2 are 'center', 'centerd', 'left', 'right',
% 'top', and 'bottom'.
%      r=AlignRect(r,screenRect,'center','top');
% For backward compatibility, also accepts 1,2,3,4,5,6. You may use the
% pre-defined constants RectLeft,RectRight,RectTop,RectBottom, but there is
% no named constant for centering and centering keeping decimal precision.
%
% rect and fixedRect can both be Mx4 or 4xM rect arrays, but either for one of
% them M must be 1, or rect and fixedRect must have the same shape. If rect
% contains multiple rects and fixedRect only one, each rect is aligned in
% fixedRect. If fixedRect contains multiple rects, but rect only one, rect
% is aligned in each of the fixedRects. If both rect and fixedRect are
% arrays with the same number of rects, each rect is aligned in the
% corresponding fixedrect. The same sequence of alignment operations will
% be executed on all rects in the array.
%
% The output will have the same orientation as rect.
%
% See also PsychRects/Contents.

% Denis Pelli 5/27/96, 7/10/96, 8/5/96, 11/8/06
% dcn 7/26/2015: vectorized, adding 'centerd' action
% dcn 7/29/2015: now handles 4xM and Mx4 inputs

if nargin<3
    error('Usage:  rect=AlignRect(rect,fixedRect,side1,[side2])');
end
if (size(rect,2)==4 && size(fixedRect,2)~=4) || (size(rect,2)~=4 && size(fixedRect,2)==4)
    % orientation of first rect(-array) is leading for orientation of
    % output. Make sure second input has same orientation
    fixedRect = fixedRect.';
end
side{1}=side1;
if nargin>3
    side{2}=side2;
end
for i=1:length(side)
    if isnumeric(side{i})
        if ~ismember(side{i},1:6)
            error('Illegal side%d value.',i);
        end
    elseif ischar(side{i})
        [ok,n]=ismember(side{i},{'left','top','right','bottom','center','centerd'});
        if ~ok
            error('Illegal side%d value.',i);
        end
        side{i}=n;
    else
        error('Illegal side%d value.',i);
    end
    if (size(rect, 1)==4) && (size(rect,2)>=1)
        % 4xM
        switch side{i}
            case 6
                rect=CenterRectd(rect,fixedRect);
            case 5
                rect=CenterRect(rect,fixedRect);
            case {1,3},
                rect=OffsetRect(rect,fixedRect(side{i},:)-rect(side{i},:),0);
            case {2,4},
                rect=OffsetRect(rect,0,fixedRect(side{i},:)-rect(side{i},:));
        end
    elseif (size(rect, 2)==4) && (size(rect,1)>=1)
        % Mx4
        switch side{i}
            case 6
                rect=CenterRectd(rect,fixedRect);
            case 5
                rect=CenterRect(rect,fixedRect);
            case {1,3},
                rect=OffsetRect(rect,fixedRect(:,side{i})-rect(:,side{i}),0);
            case {2,4},
                rect=OffsetRect(rect,0,fixedRect(:,side{i})-rect(:,side{i}));
        end
    else
        % Something weird and unknown:
        error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
    end
end
