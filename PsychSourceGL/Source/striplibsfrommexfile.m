function striplibsfrommexfile(filename, testrun)
% striplibsfrommexfile(filename [, testrun=0]);
%
% Remove a certain set of dynamic library dependencies from MEX/OCT file
% given via 'filename'. Optional 'testrun' if set to non-zero value will
% perform a dry-run.
%
% The stripped libraries (see source of this file) are library dependencies
% added by default by Octave's build process for OCT/MEX files. They make
% sense for mathematical routines, but are totally redundant for PTB files.
%
% They also do harm, as they are versioned for some odd reason, so trying
% to load the mex files on Linux distributions with different versions of
% those libraries installed will fail.
%
% -> This script loads the binary code, searches for occurence of the
% superfluous library name strings and null's them out to remove the
% dependency. Then it rewrites the stripped mex binary.
%

if nargin < 2
    testrun = 0;
end

if exist(filename, 'file')
    % Read original image:
    [fd, msg] = fopen(filename, 'r+');
    if fd == -1
        error('ERROR! MEX file %s could not be read [%s]! Can''t strip it!\n', filename, msg);
    end
    image = uint8(fread(fd));
 
    % Strip all following libraries from image:
    image = stripLibrary(image, 'libreadline');
    image = stripLibrary(image, 'libncurses');
    image = stripLibrary(image, 'libfft');
    image = stripLibrary(image, 'liblapack');
    image = stripLibrary(image, 'libblas');
    image = stripLibrary(image, 'libhdf');
    image = stripLibrary(image, 'libgfortran');
    image = stripLibrary(image, 'libz');
    image = stripLibrary(image, 'libcruft');
    image = stripLibrary(image, 'liboctave');

    % Write stripped image:
    frewind(fd);
    if ~testrun
        fwrite(fd, image);
    end
    
    if (fclose(fd) ~= 0)
        error('ERROR! MEX file %s could not be closed [%s]! Can''t strip it!\n', filename, ferror(fd));
    end

    fprintf('MEX file %s stripped from superfluous libraries.\n', filename);
else
    fprintf('WARNING! MEX file %s does not exist! Can''t strip it!\n', filename);
end

return

function image = stripLibrary(image, library)

    dodebug = 0;

    % Find first character of library in image:
    pStart = strfind(char(image'), library);
    if ~isempty(pStart)
        % Iterate over all occurences, kill each of them:
        for cStart = pStart
            % Kill occurence at cStart:
            
            if dodebug
                disp(char(image(cStart:cStart + 20)))
            end

            % Start at the start:
            p = cStart;

            % Strip-Loop: Strip until zero-byte:
            while (image(p) ~= 0)
                % Null-out this character:
                image(p) = 0;
                % Next one...
                p = p + 1;
            end

            if dodebug
                disp(char(image(cStart:cStart + 20)))
            end
        end
    end
    
    if ~isa(image, 'uint8')
        error('After strip op, image is no longer of uint8 class as required! Stripping failed!!');
    end
    
return
