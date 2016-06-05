function osxsetoctaverpath(mexfname, mexpath)
% osxsetoctaverpath(mexfname [, mexpath])
%
% Change the @rpath library search path for the octave
% runtime libraries inside the given mex file.
%
% As the mex files for Octave-4 are still built on Octave-3,
% we change from absolute path to rpath and from the version 2
% libraries to version 3 libraries.
%
% E.g.,
%
% osxsetoctaverpath('Screen'); would rewrite Screen.mex
% to use the @rpath settings stored in this function and
% work on Octave-4 instead of Octave-3.
%

    if ~IsOSX(1) || ~IsOctave
        error('osxsetoctaverpath only works with a 64-Bit version of Octave for OSX!');
    end

    % If no mex filename given, iterate over 'mexpath' - or the default install
    % location of mex files - and apply the rpath editing to each mex file there:
    if nargin < 1 || isempty(mexfname)
        if nargin < 2 || isempty(mexpath)
            mexpath = [PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/'];
        end

        d = dir (mexpath);
        for j = 1:length(d)
            if ~d(j).isdir
                [a, mexfname, extension] = fileparts(d(j).name);
                if ~isempty(strfind(extension, mexext))
                    osxsetoctaverpath(mexfname, mexpath);
                end
            end
        end
        return;
    end

    % Set default path for finding the mex file to process, if omitted:
    if nargin < 2 || isempty(mexpath)
        mexpath = '../Projects/MacOSX/build/';
    end

    % Build full path to file:
    mexfname = [mexpath mexfname '.' mexext];

    % Replace absolute path to liboctinterp.2.dylib with @rpath:
    system(['install_name_tool -change ' octave_config_info.octlibdir '/liboctinterp.2.dylib @rpath/liboctinterp.3.dylib ' mexfname]);

    % Replace absolute path to liboctave.2.dylib with @rpath:
    system(['install_name_tool -change ' octave_config_info.octlibdir '/liboctave.2.dylib @rpath/liboctave.3.dylib ' mexfname]);

    % Add one single rpath: @loader_path. This is the path to our folder where our
    % mex file is stored. If we place copies of liboctave.3.dylib and liboctinterp.3.dylib
    % there, then the linker will find them. In absence, the linker will also search the
    % users $HOME/lib/ directory as a possible fallback:
    lpaths = { '@loader_path' };

    % Add all paths in lpaths as potential search paths for the octave
    % library directories, ie., as settings for @rpath:
    for i = 1:length(lpaths)
        system(['install_name_tool -add_rpath ' lpaths{i} ' ' mexfname]);
        fprintf('Added Octave-4 @rpath %s to mex file %s ...\n', lpaths{i}, mexfname);
    end

return;
