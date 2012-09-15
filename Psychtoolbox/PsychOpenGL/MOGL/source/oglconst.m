function oglconst(glheaderpath, aglheaderpath)

% OGLCONST  Collect GL, GLU, and AGL constants from C header files, and
%           store them in oglconst.mat
%
% usage:  oglconst

% 09-Dec-2005 -- created (RFM)
% 23-Jan-2005 -- constants saved in both struct and OpenGL style (RFM)
% 05-Mar-2006 -- ability to spec. system header file path added (MK)
% 30-Aug-2012 -- Also parse own glext.h file for up-to-date definitions (MK)
% 30-Aug-2012 -- Make more robust: Handle 64-Bit integer defines and their
%                quirks, as well ad #define GL_THIS GL_THAT "recursive"
%                definitions (MK).

if IsWin
    error('Parsing of GL header files on Windows not yet supported.');
end;

% Alternate path to header files specified?
if nargin < 1
    if IsOSX
        glheaderpath = '/System/Library/Frameworks/OpenGL.framework/Headers';
    end;

    if IsLinux
        glheaderpath = '/usr/include/GL';
    end;
end;

if nargin < 2
    if IsOSX
        aglheaderpath = '/System/Library/Frameworks/AGL.framework/Headers';
    end;
end;

if IsOSX
    fprintf('Parsing OpenGL, GLU and AGL header files in %s and %s ...\n',glheaderpath, aglheaderpath);
else
    fprintf('Parsing OpenGL and GLU header files in %s ...\n',glheaderpath);
end;

% get constants from the OpenGL header files.  the 'parsefile' routine
% defines the constants in the calling workspace, as variables with the
% same names as the #defined OpenGL constants, e.g., GL_COLOR_BUFFER_BIT.
% the return argument contains all the constants as fields of a structure,
% e.g., GL.COLOR_BUFFER_BIT.  (note that the first underscore has been
% changed to a period.)  the first style is more like the C convention,
% and the second style prevents the workspace from being swamped with
% hundreds of global variables.  later on, we can load whichever style
% we want from the file where they're all saved.
GL= parsefile(sprintf('%s/gl.h', glheaderpath), 'GL_');
GL= parsefile(sprintf('%s/glext.h', glheaderpath), 'GL_', GL);

% Also parse a glext.h file if it is located in our headers/ subdirectory
% under glext.h
if exist('./headers/glext.h', 'file')
    fprintf('Parsing also our private glext.h file for more up-to-date constant definitions.\n');
    GL= parsefile('./headers/glext.h', 'GL_', GL);
end

GLU=parsefile(sprintf('%s/glu.h', glheaderpath),'GLU_'); %#ok<*NASGU>
if IsOSX
    AGL=parsefile(sprintf('%s/agl.h', aglheaderpath),'AGL_');
end;

fname='oglconst.mat';

% save OpenGL-style constants
if IsOSX
    save(fname,'GL_*','GLU_*','AGL_*', '-V6');
    % save structure-style constants to same file
    save(fname,'GL','GLU','AGL','-append', '-V6');
end;

if IsLinux
    save(fname,'GL_*','GLU_*','-V6');
    % save structure-style constants to same file
    save(fname,'GL','GLU','-append','-V6');
end;

% put a copy into the 'core' directory
copyfile(fname,'../core');

return


% function to parse header files
function S = parsefile( fname, prefix, origin )

% initialize return argument
if nargin < 3
    S=[];
else
    S=origin;
end;

% check size of prefix (GL, GLU, or AGL)
nprefix=length(prefix);

% open input file
fid=fopen(fname,'r');

% step through lines of input file
while ~feof(fid),

    % read a line and parse it as '#define SYMBOL VALUE'
    codeline=fgets(fid);
    r=regexp(codeline,'^\s*#define\s+(?<symbol>\S+)\s*(?<value>\S*)','names');

    % if it's not a #define statement, then skip it
    if isempty(r),
        continue
    end

    % if the symbol doesn't have the required prefix, (GL_, GLU_, or AGL_),
	% then skip it
    if ~strncmp(r.symbol,prefix,nprefix),
        continue
    end

    % remove prefix from symbol
	fieldname=r.symbol((nprefix+1):end);
    
    % if remainder of symbol name begins with a digit, then add 'N' to make
    % it a valid field name
    if ismember(fieldname(1),'0123456789'),
        fieldname=[ 'N' fieldname ]; %#ok<*AGROW>
    end

    % convert value to a numeric value
    if ~isempty(r.value),
        try
            if strncmp(r.value,'0x',2),
                % convert hex value:
                
                % First strip trailing u or l characters, which would
                % define (u)nsigned (l)ong values etc., as hex2dec can't
                % handle such trailing specifiers:
                inpstring = r.value(3:end);
                inpstring = inpstring(inpstring ~= 'u');
                inpstring = inpstring(inpstring ~= 'l');
                
                % Convert filtered string to decimal:
                nvalue=hex2dec(inpstring);
            elseif strncmp(r.value,prefix,length(prefix)),
                % Constant defined by name of a previously defined
                % constant, e.g., #define GL_THIS GL_THAT.
                %
                % Try to find value of predecessor constant:
                substitutionField = r.value(length(prefix)+1:end);
                if isfield(S, substitutionField)
                    nvalue = getfield(S, substitutionField); %#ok<*GFLD>
                    fprintf('\n-> Substituted %s with %s value of %x \n', fieldname, substitutionField, nvalue);
                else
                    fprintf('\n-> FAILED Substitution of %s with %s due to missing preceeding definition!\n', fieldname, substitutionField);
                    nvalue = [];
                end
            else
                % convert decimal value
                nvalue=str2num(r.value); %#ok<*ST2NM>
            end
        catch %#ok<*CTCH>
            nvalue = [];
        end
        % assign value of zero if conversion failed, e.g., if symbol was
		% defined using another #define, as in #define GL_THIS GL_THAT.
		% if anything important is #defined this way, we'll have to fix
		% up this part.
        if isempty(nvalue),
            warning('error converting numeric value from line:  %s',codeline); %#ok<*WNTAG>
            nvalue=0;
        end
	else
		% assign zero if symbol has no value, e.g., #define GL_I_WAS_HERE
        nvalue=0;
    end

    % add numeric value to struct
    fprintf(1,'%s\t%-20s\t%s\t%.0f\n',prefix,fieldname,r.value,nvalue);
    S=setfield(S,fieldname,nvalue); %#ok<*SFLD>
    
    % define OpenGL-style variable in calling workspace
	evalin('caller',sprintf('%s=%d;',r.symbol,nvalue));

end

% close file
fclose(fid);

return
