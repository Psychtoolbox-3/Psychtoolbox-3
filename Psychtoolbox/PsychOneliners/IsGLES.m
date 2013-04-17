function rc = IsGLES
% Returns 1 if active rendering api is OpenGL-ES,
% 0 otherwise.

val = getenv('PSYCH_USE_GFX_BACKEND');
rc = (length(val) >= 4) && strcmp(val(1:4), 'gles');
