function retval = moglsingle(arg)
% --- OBSOLETE --- OBSOLETE --- OBSOLETE --- OBSOLETE ---
% PSYCHTOOLBOX SPECIFIC single() implementation:
%
% retval = moglsingle(arg) -- convert into
% a single precision floating point number.
%
% This routine takes an 'arg' of arbitrary
% numeric class and converts it into an equivalent
% object of single precision floating point format.
%
% If a builtin single() function is available,
% as on Matlab and Octave 3.2+, it calls the builtin
% single() function.
%
% Otherwise (Octave) it would call our own special
% implementation.
%
% This is no longer needed as of Octave 3.2.0, but we leave the function
% here as many internal and external code relies on its presence.
%
% For Octave pre 3.2, this applied:
%
% This is a hack needed to make OpenGL (MOGL) work
% on GNU/Octave, despite Octave's lack of a single
% precision data type. Returned values are not
% useful or accessible by Octave. They only have
% meaning for the Psychtoolbox OpenGL functions.
%

% ---protected---
retval = single(arg);
