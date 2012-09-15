% Psychtoolbox:PsychDemos:MovieDemos
%
% Demos demonstrating use of the Quicktime and GStreamer movie playback
% functions. Quicktime is the default engine on 32-Bit OS/X. It is only
% available on 32-Bit OS/X. GStreamer is the default engine on Linux,
% Windows, and 64-Bit OS/X.
%
% All demos work with the GStreamer engine. Most demos also work with
% Quicktime, although Quicktime is deprecated and will no longer be
% improved any further.
%
% Read the help for the corresponding Screen subfunctions for infos.
%
% SimpleMovieDemo              -- Most simplistic demo on how to play a movie.
% PlayMoviesDemo               -- Show simple playback of one movie with sound at a time.
% PlayDualMoviesDemo           -- Same as PlayMoviesDemo, but play two movies in parallel.
% PlayInterlacedMovieDemo      -- Simple demo on how to use the GLSL realtime video deinterlacer.
% PlayMoviesWithoutGapDemo1    -- Play one movie while opening another one for gap-less playback.
% PlayMoviesWithoutGapDemo2    -- GStreamer only: Optimized version of gap-less playback.
%                                 movie in the background to reduce gaps between movies.
% DetectionRTInVideoDemo       -- How to collect reaction times in response to detection of some
%                                 event in a presented movie file. Takes care to get timing right.
% LoadMovieIntoTexturesDemo    -- Quickly load a movie into a stack of textures for quick playback
%                                 with arbitrary speed and order.
%
