function autocode(overwrite, glheaderpath, openal)

% AUTOCODE.M  Generate MATLAB-OpenGL or OpenAL interface code from OpenGL or OpenAL header files
%
% Usage: autocode(overwrite, glheaderpath, openal).
% glheaderpath is the path to the OpenGL gl.h and glu.h header files.
% It defaults to /System/Library/Frameworks/OpenGL.framework/Headers/
%
% If overwrite is == 1, then existing M-Files are overwritten, else (and
% this is the default) existing files are not touched.
%
% If openal is == 1, then OpenAL code is generated, instead of OpenGL code.
%
% 18-Dec-05 -- created (RFM)
% 05-Mar-06 -- added option to spec. glheaderpath (MK)
% 30-May-06 -- added option 'overwrite' to not overwrite M-Files. (MK)
%              added additional parsing code for headers/glext_edit.h
% 06-Feb-07 -- added support for OpenAL (MK)
% 24-Mar-11 -- Silence mlint warnings. (MK)
% 01-Apr-12 -- Adapt to parsing of current glext_edit.h from OpenGL registry. (MK)

clc;

if nargin < 1 || isempty(overwrite)
    overwrite = 0;
end;
overwrite %#ok<NOPRT>

if nargin < 3
	openal = 0;
end

% Alternate path to header files specified?
if nargin < 2 || isempty(glheaderpath)
    if IsOSX
        glheaderpath = '/System/Library/Frameworks/OpenGL.framework/Headers';
		if openal
			glheaderpath = '/System/Library/Frameworks/OpenAL.framework/Headers';
		end
    end;
    if IsLinux
        glheaderpath = '/usr/include/GL';
		if openal
			glheaderpath = '/usr/include/AL';
		end
    end;
    if IsWin
        error('autocode.m not yet supported on M$-Windows.');
    end;
end;

if openal
	fprintf('Parsing OpenAL and ALC header files in %s ...\n',glheaderpath);
else
	fprintf('Parsing OpenGL and GLU header files in %s ...\n',glheaderpath);
end

% delete unprotected wrapper files
% ( protected = contains string "---protected---" somewhere in the file )
if overwrite > 0
    ! find ../wrap/*.m -not -exec grep -q -e ---protected--- {} \; -print0 | xargs -0 -n 1 rm
end;

% make file with list of OpenGL functions
tmplistfile='/tmp/mogl_listfile.txt';

if openal
	unix(sprintf('grep al[A-Z]   %s/al.h        | grep -v ProcPtr | grep -v \\#define | sed -E ''s/^[[:space:]]*extern[[:space:]]*//'' | sed -E ''s/;[[:space:]]*$//'' >  %s',glheaderpath, tmplistfile));
	unix(sprintf('grep alu[A-Z]  %s/alc.h       | grep -v ProcPtr | grep -v \\#define | sed -E ''s/^[[:space:]]*extern[[:space:]]*//'' | sed -E ''s/;[[:space:]]*$//'' >> %s',glheaderpath, tmplistfile));

	% initialize and open C file
	cfile='al_auto.c';
	if unix(sprintf('cat al_auto_init.c | sed ''s/DATE/%s/'' > %s',datestr(now,1),cfile))~=0,
		error(sprintf('unable to initialize ''%s''',cfile)); %#ok<SPERR>
	end
else
	unix(sprintf('grep gl[A-Z]   %s/gl.h        | grep -v ProcPtr | grep -v \\#define | sed -E ''s/^[[:space:]]*extern[[:space:]]*//'' | sed -E ''s/;[[:space:]]*$//'' >  %s',glheaderpath, tmplistfile));
	unix(sprintf('grep glu[A-Z]  %s/glu.h       | grep -v ProcPtr | grep -v \\#define | sed -E ''s/^[[:space:]]*extern[[:space:]]*//'' | sed -E ''s/;[[:space:]]*$//'' >> %s',glheaderpath, tmplistfile));
	unix(sprintf('grep glut[A-Z] headers/glut_edit.h | grep -v ProcPtr | grep -v \\#define | sed -E ''s/^[[:space:]]*extern[[:space:]]*//'' | sed -E ''s/;[[:space:]]*$//'' >> %s',tmplistfile));
	unix(sprintf('grep gl[A-Z]   headers/glext_edit.h | grep -v ProcPtr | grep -v \\#define | sed -E ''s/^[[:space:]]*GLAPI[[:space:]]*//''  | sed -E ''s/[[:space:]]*APIENTRY*//'' | sed -E ''s/;[[:space:]]*$//'' >>  %s',tmplistfile));
    
	% initialize and open C file
	cfile='gl_auto.c';
	if unix(sprintf('cat gl_auto_init.c | sed ''s/DATE/%s/'' > %s',datestr(now,1),cfile))~=0,
		error(sprintf('unable to initialize ''%s''',cfile)); %#ok<SPERR>
	end
end

listfid=fopen(tmplistfile,'r');
cfid=fopen(cfile,'a');

% initialize list of entries in command map
cmdmap={};

% step through list of OpenGL functions
while(~feof(listfid)),

    % parse the next function declaration
    funcp=cparse(fgetl(listfid), openal);

    % skip functions in autono.txt
    if isempty(funcp) || filecontains('autono.txt',funcp.fname),
        continue
    end

    % skip functions that use double indirection in an input argument
    if doubleindirect(funcp),
        % add function name to autono.txt
        unix(sprintf('cat "%s  %% uses double indirection (added to autono.txt by autocode.m) >> autono.txt',funcp.fname));
        continue
    end

    % function meets criteria, so show function name
    fprintf(1,'Found: %s ...',funcp.fname);

    % extract information needed to construct M-file wrapper and
    % C interface function (input arguments, return arguments, etc.)
    [M,C]=mcinfo(funcp);

    if (overwrite > 0) || (exist([ funcp.fname '.m'], 'file')~= 2)
        % (re)create M-file wrapper
        fprintf(1, 'creating M-File wrapper.\n');
        mwrite(funcp,M,openal);
    else
        fprintf(1, '\n');
    end;
    
    % create C interface function
	cwrite(cfid,funcp,C);
	
	% add line to command map (to be written to C file later)
	cmdmap{end+1}=sprintf('{ %-35s%-35s }',[ '"' funcp.fname '",' ],C.interfacename); %#ok<AGROW>
	
end

% sort and write command map
cmdmap=sort(cmdmap);
fprintf(cfid,'int gl_auto_map_count=%d;\n',numel(cmdmap));
fprintf(cfid,'cmdhandler gl_auto_map[] = {\n');
fprintf(cfid,'%s,\n',cmdmap{1:end-1});
fprintf(cfid,'%s};\n',cmdmap{end});

% close files
fclose(cfid);
fclose(listfid);

% clean up
delete('/tmp/mogl_*');
delete('/tmp/moal_*');

return
