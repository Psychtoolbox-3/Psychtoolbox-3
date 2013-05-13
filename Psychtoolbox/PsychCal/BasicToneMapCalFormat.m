function outputXYZCalFormat = BasicToneMapCalFormat(inputXYZCalFormat, maxLum)
% outputXYZCalFormat = BasicToneMapCalFormat(inputXYZCalFormat, maxLum)
%
% Simple tone mapping.  Leaves any pixel with luminance below maxLum alone.
% For pixels whose luminance exceeds maxLum, scale XYZ down multiplicatively so
% that luminance is maxLum.
%
% 10/1/09 bjh, dhb     Created it.
% 10/4/09 dhb          Debug and make it work right.

% Find offending pixels
index = find(inputXYZCalFormat(2,:) > maxLum);

% If any pixel exceeds maxLum, scale it by 1/Y.  Uses
% MATLAB's indexing trick of repeating an index to replicate 
% values.
outputXYZCalFormat = inputXYZCalFormat;
if (~isempty(index))
   outputXYZCalFormat(:, index) = maxLum*(inputXYZCalFormat(:, index)./inputXYZCalFormat([2 2 2]', index)); 
end

end
