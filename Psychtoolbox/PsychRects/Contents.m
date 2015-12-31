% Psychtoolbox:PsychRects.
%
% help Psychtoolbox % For an overview, triple-click me & hit enter.
%
%     AdjoinRect          - Move one rect next to another.
%     AlignRect           - Align one rect over another.
%     ArrangeRects        - Arrange array of rects in a pleasing way.
%     CenterRect          - Center the first rect in the second.
%     CenterRectArrayOnPoint - Centers a group of rectangles 'rects' on a given point.
%     CenterRectOnPoint   - Center rect around an x,y coordinate, round to
%                           integral positions. Vectorized.
%     CenterRectOnPointd  - Center rect around an x,y coordinate. Vectorized.
%     ClipRect            - Return intersection of two rects.
%     GrowRect            - Grow a rect by a certain number of pixels.
%     InsetRect           - Shrink/expand rect by additive insets.
%     IsEmptyRect         - Returns 1 if empty, returns 0 otherwise.
%     IsInRect            - Is the point inside the rect?
%     OffsetRect          - Shift rect vertically and horizontally. Vectorized.
%     RectBottom          - Return index of yBottom entry of a rect.
%     RectCenter          - Returns the integer x,y coordinates of center. Vectorized.
%     RectCenterd         - Returns the exact x,y coordinates of center. Vectorized.
%     RectOfMatrix        - Accept a matrix and return its bounding rect.
%     RectHeight          - Return height of a rect.
%     RectLeft            - Return index of xLeft entry of a rect.
%     RectRight           - Return index of xRight entry of a rect.
%     RectTop             - Return index of yTop entry of a rect.
%     RectWidth           - Return width of a rect.
%     RectSize            - Return the width and the height of a rect.
%     ScaleRect           - Scale a rect by multiplicative factors.
%     SetRect             - Create a rect.
%     SizeOfRect          - Accept a Psychtoolbox rect and return the
%                           corresponding image's MATLAB array size.
%     UnionRect           - Smallest rect containing two given rects.
%
% In general, inputs for most functions can be rects or rect-arrays. Rect
% arrays are 4xM or Mx4 matrices in which each row is a rect. The function
% help will clarify how rect arrays are treated in the function.
% NB: 4x4 inputs will be treated as 4xM input, which means it is not
% advisable to use Mx4 rect arrays unless you are sure that M will never be
% 4.

% Copyright (c) 1997-2015 by David Brainard & Denis Pelli and the PTB
% developers.
