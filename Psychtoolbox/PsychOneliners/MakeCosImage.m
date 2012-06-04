function [image] = MakeCosImage(freqi,freqj,nRowPixels,nColPixels)
% [image] = MakeCosImage(freqi,freqj,nRowPixels,[nColPixels])
%
% Computes a two-dimensional cosine function image.
%
% The image has dimensions nRowPixels by nColPixels.
% If nColPixels is omitted, a square image is returned.
%
% 8/15/94		dhb		Both row and column dimensions used if passed.

% Set column pixels for square image if it wasn't passed.
if (nargin <= 3)
	nColPixels = nRowPixels;
end

x = 1:nColPixels;
y = 1:nRowPixels;
usefreqi = 2*pi*freqi/nRowPixels;
usefreqj = 2*pi*freqj/nColPixels;
cosx = cos(usefreqj*x);
cosy = cos(usefreqi*y);
image = cosy'*cosx;




