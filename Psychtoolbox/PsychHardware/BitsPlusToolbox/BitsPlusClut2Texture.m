function texturePtr = BitsPlusClut2Texture(windowPtr, clut)
% texturePtr = BitsPlusClut2Texture(windowPtr, clut)
%
%   Generates a PsychToolbox texture containing the CLUT + magic code
%   required to set the clut in Bits++ mode.
%
%   'clut' should be a 256x3 matrix consisting of values in the range
%   specified by Screen('ColorRange'), which is by default [0, 255].

if nargin ~= 2
    error('Usage: texturePtr = BitsPlusClut2Texture(windowPtr, clut)');
end

% Convert the clut into the current color range.
colorRange = Screen('ColorRange', windowPtr);
clut = clut / colorRange .* (2^16 - 1);

% Convert the clut into the special one for Bits++.
newClutRow = BitsPlusEncodeClutRow(clut);

% Generate the texture holding the Bits++ clut.
texturePtr = Screen('MakeTexture', windowPtr, uint8(newClutRow));
