function rc = IsGLES1
% Returns 1 if active rendering api is OpenGL-ES1.x,
% 0 otherwise.

rc = strcmp(getenv('PSYCH_USE_GFX_BACKEND'), 'gles1');
