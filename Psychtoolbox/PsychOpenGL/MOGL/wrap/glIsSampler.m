function r = glIsSampler( sampler )

% glIsSampler  Interface to OpenGL function glIsSampler
%
% usage:  r = glIsSampler( sampler )
%
% C function:  GLboolean glIsSampler(GLuint sampler)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

r = moglcore( 'glIsSampler', sampler );

return
