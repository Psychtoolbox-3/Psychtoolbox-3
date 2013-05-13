function CloseOpenAL
% CloseOpenAL - Shutdown OpenAL based audio subsystem.
%
% Use this at the end of your script in the cleanup- or error-handling
% routine to stop and shutdown an OpenAL audio device which was previously
% initialized via InitializeMatlabOpenAL.

% History:
% 7.2.2007 Written (MK).

moalcore('DEBUGLEVEL', 1);
moalcore('almClose')
return
