function rect = RectOfMatrix(mat)
% rect = RectOfMatrix(mat)
%
% Accept an image as a matrix and return a Psychtoolbox rect specifying 
% the bounds.
%
% See also PsychRects/Contents, SetRect.

% HISTORY
% 6/28/04   awi Wrote it.
% 7/13/04   awi Improved documentation.  
% 1/29/05   dgp Renamed from RectFromMatrix to RectOfMatrix. 
% 1/29/05   dgp Timed it. RectOfMatrix is quick, taking about 100 us, 
%               independent of matrix size.
% 11/24/05  mk  Bugfix. Didn't work for multi-layer (3D) matrices.

[m,n,~]=size(mat);
rect=[0 0 n m];
