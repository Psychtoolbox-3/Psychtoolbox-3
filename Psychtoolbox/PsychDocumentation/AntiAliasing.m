% AntiAliasing -- Psychtoolbox facilities for Anti-Aliasing.
%
% According to http://www.wikipedia.org , anti-aliasing is the
% technique of minimizing aliasing (jagged or blocky patterns)
% when representing a high-resolution signal at a lower resolution.
%
% If you want to know in detail about the reason for aliasing artifacts
% in digital computer images and the theoretical approaches for removing
% or reducing it, as a starter, search and read on Wikipedia about "Aliasing"
% and "Anti-Aliasing", or read a book or website about computer graphics
% or digital image processing.
%
% Psychtoolbox allows for multiple approaches to anti-alias your stimuli.
% The methods all represent some tradeoff between ease of use, generality
% and amount of brute-force applied to implement them.
%
% Multisampling / Supersampling:
% 
% The easiest (for you) and at the same time most general approach is the
% brute-force approach, called "Multisampling" or "Supersampling", also
% known as "Full scene anti-aliasing" or FSAA. You can enable multisampling
% for an onscreen window by simply setting the optional argument 'multisample'
% of the Screen('OpenWindow', ...) - subfunction to a value greater than zero.
% Psychtoolbox will check if your graphics hardware supports multisampling and
% enable it for the new onscreen window, if possible. In a nutshell, the graphics
% hardware will then perform all drawing operations (including all images, all 2D shapes
% like points, discs, lines, rectangles and polygons, text and all 3D OpenGL rendering)
% at a much higher resolution than your display resolution, computing multiple color
% samples for each output pixel. When calling Screen('Flip'), the hardware will
% downsample this high resolution image by use of a proper lowpass filter, combining
% the multiple color values of the multiple samples per pixel into a single color
% value for the final output pixel. A typical algorithm for this would be a weighted
% average of the color values of the multiple samples, but depending on your gfx-
% hardware, the algorithm may use more clever schemes. The value of the 'multisample'
% parameter roughly specifies, how many samples should be used for each final pixel.
% A higher value means a higher quality, but at the same time higher requirements in
% terms of video memory usage and processing time.
%
% Psychtoolbox tries to select a value as close as possible to the requested value: It
% rounds up to the closest value supported by your hardware (e.g., if you request 6
% samples, but the hardware only supports either 4 or 8 samples, PTB will select 8
% samples). If you request a value that is greater than the highest supported value,
% PTB will use the maximum value supported by your system.
%
% The algorithms used for multisampling differ between each model of graphics hardware:
% NVidia's gfx-cards do it differently than ATI's cards and each model from NVidia or
% ATI also does it differently from its predecessors. For a general overview of how
% multisampling and supersampling works, have a look at Wikipedias articles on
% "Supersampling". For specific gfx-cards, the hardware vendors provide white papers
% and specifications on how their devices operate at specific settings, e.g., one
% presentation of NVidia is http://developer.nvidia.com/object/gdc_ogl_multisample.html
% A more detailed article is http://www.3dcenter.org/artikel/geforcefx_aa_modi/index_english.php
%
% As stated above, this brute-force approach of anti-aliasing requires minimal
% work for you (setting one optional parameter) and it works always, but it is
% expensive. For a multisample value of n, your hardware will need roughly n+1 times
% the amount of video memory and it will need more than n times the amount of computation
% time for drawing and showing your stimulus. The maximum value that you can set for n
% therefore not only depends on how modern your graphics-card is, but is also limited
% by the amount of VRAM installed on your gfx-card, the screen resolution you've selected,
% the complexity of your stimulus and the stimulus update rate you need, so choose a
% reasonable tradeoff for your purpose.
%
% Some sort of multisampling/supersampling should be supported by basically any piece
% of graphics hardware for Intel PC's under M$-Windows or GNU/Linux sold by NVidia or ATI 
% after beginning of 2001. On MacOS-X, all supported gfx-adapters, except the ATI Rage 128
% and ATI Rage 128 Pro, should allow for at least 4-way sampling.
%
% Typical values on todays hardware: GeForce-3 series: 4-way sampling,
% GeForce 7800 series: 16 way multisampling.
%
% To be continued...
