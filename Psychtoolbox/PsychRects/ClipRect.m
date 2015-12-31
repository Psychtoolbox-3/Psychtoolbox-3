function newRect = ClipRect(a,b)
% newRect = ClipRect(a,b)
%
% Returns the rect that is the intersection of the two rects a and b.
% Returns an empty rect [0 0 0 0] if the two rects are apart.
%
% a and b can each be single rects or Mx4 or 4xM rect-arrays. If both are
% rect-arrays, they must contain the same number of rects.
%
% The output will have the same orientation as a.
%
% See also PsychRects/Contents.

% 7/5/96  dgp  Wrote it.
% 7/6/96  dgp  Return empty matrix [] only if apart.
% 7/10/96 dgp  Return empty rect [0 0 0 0] if apart.
% 7/26/15 dcn  Vectorized.
% 7/29/15 dcn  Now handles 4xM and Mx4 inputs

if nargin~=2
    error('Usage:  rect=ClipRect(a,b)');
end
if (size(a,2)==4 && size(b,2)~=4) || (size(a,2)~=4 && size(b,2)==4)
    % orientation of first rect(-array) is leading for orientation of
    % output. Make sure second input has same orientation
    b = b.';
end

if (size(a, 1)==4) && (size(a,2)>=1)
    % 4xM
    if size(b,2) > size(a,2)
        newRect=repmat(a,1,size(b,2)./size(a,2));
    else
        newRect = a;
    end
    
    newRect(2,:) = max(a(2,:) ,b(2,:));
    newRect(4,:) = min(a(4,:) ,b(4,:));
    newRect(1,:) = max(a(1,:) ,b(1,:));
    newRect(3,:) = min(a(3,:) ,b(3,:));
    
    qZero = IsEmptyRect(newRect);
    newRect(:,qZero) = 0;
elseif (size(a, 2)==4) && (size(a,1)>=1)
    % Mx4
    if size(b,1) > size(a,1)
        newRect=repmat(a,size(b,1)./size(a,1),1);
    else
        newRect = a;
    end
    
    newRect(:,2) = max(a(:,2) ,b(:,2));
    newRect(:,4) = min(a(:,4) ,b(:,4));
    newRect(:,1) = max(a(:,1) ,b(:,1));
    newRect(:,3) = min(a(:,3) ,b(:,3));
    
    qZero = IsEmptyRect(newRect);
    newRect(qZero,:) = 0;
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end
