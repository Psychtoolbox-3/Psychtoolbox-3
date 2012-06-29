function [output] = Trunc(input)
% [output] = Trunc(input)
% Truncate to range [0-1].

% 6/13/12    dn     Use boolean indices, scalars expand to all needed outputs

output = input;

qIndex = input < 0;
if any(qIndex(:))
    output(qIndex) = 0;
end

qIndex = input > 1;
if any(qIndex(:))
    output(qIndex) = 1;
end
