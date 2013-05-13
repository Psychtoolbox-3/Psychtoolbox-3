function [image] = MakeSineImage(freqi,freqj,nRowPixels,nColPixels)
% [image] = MakeSineImage(freqi,freqj,nRowPixels,[nColPixels])
%
% Computes a two-dimensional sine function image.
%
% The image has dimensions nRowPixels by nColPixels.
% If nColPixels is omitted, a square image is returned.
%
% 8/15/94		dhb		Both row and column dimensions used if passed.
%				dhb		Changed zero frequency convention.
% 6/20/98       dhb, mw Fixed error in zero handling case.

% Set column pixels for square image if it wasn't passed.
if (nargin <= 3)
	nColPixels = nRowPixels;
end

x = 1:nColPixels;
y = 1:nRowPixels;
usefreqi = 2*pi*freqi/nRowPixels;
usefreqj = 2*pi*freqj/nColPixels;

% Handle zero frequency case
if (usefreqj == 0 && usefreqi ~= 0)
	sinx = ones(size(x));
else
	sinx = sin(usefreqj*x);
end

% Handle zero frequency case
if (usefreqi == 0 && usefreqj ~= 0)
	siny = ones(size(y));
else
	siny = sin(usefreqi*y);
end

% Build composite image
image = siny'*sinx;
