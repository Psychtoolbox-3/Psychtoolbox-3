function p = cparse( str , openal)

% CPARSE  Parse a C language function declaration
% 
% p = cparse( str )

% 22-Dec-2005 -- created (RFM)
% 28-Aug-2012 -- Modified; Handle datatypes GLint64, GLuint64 and GLsync. (MK)
% 28-Aug-2012 -- Modified; Make more robust against unparseable functions. (MK)

if nargin<2
    openal = 0;
end

%str = [];

if nargin<1 || isempty(str),
    if openal
        str='  const AL_API void AL_APIENTRY alFunction (  int arg1 , const char** arg2 ) ; ';
    else 
        str='  const  void * *glFunction (  int arg1 , const char** arg2 ) ; ';
    end
end

% delete extra spaces
str=sed(str,'s/^[[:space:]]*//','s/[[:space:]]*$//','s/[[:space:]]+/ /g');

% record full declaration
p.full=str;

% make list of data types
ctypes='void|char|unsigned char|signed char|short|unsigned short|signed short|short int|unsigned short int|signed short int|int|unsigned int|signed int|long|unsigned long|signed long|long int|unsigned long int|signed long int|float|double';
ogltypes='GLenum|GLboolean|GLbitfield|GLbyte|GLshort|GLint64|GLint|GLsizei|GLubyte|GLushort|GLuint64|GLuint|GLfloat|GLclampf|GLdouble|GLclampd|GLvoid|GLintptr|GLsizeiptr|GLchar|GLUnurbs|GLUquadric|GLUtesselator|GLhandleARB|GLsync';
oaltypes='ALenum|ALboolean|ALbitfield|ALbyte|ALshort|ALint|ALsizei|ALubyte|ALushort|ALuint|ALfloat|ALclampf|ALdouble|ALclampd|ALvoid|ALintptr|ALsizeiptr|ALchar';
types=[ ctypes '|' ogltypes '|' oaltypes];

% parse full declaration
if openal
    % OpenAL parser expression: ?(?<stars>[\\* ]*)
    r=regexp(str,sprintf('(?<const>const)? .* ?(?<basetype>%s) .*  ?(?<fname>\\w+) ?\\( ?(?<argin>.*) ?\\).*',types),'names');
    if ~isempty(r)
        r.stars = '';
    end
else
    % OpenGL parser expression:
    r=regexp(str,sprintf('(?<const>const)? ?(?<basetype>%s) ?(?<stars>[\\* ]*) ?(?<fname>\\w+) ?\\( ?(?<argin>.*) ?\\).*',types),'names');
end

if isempty(r) || length(r) < 1   
    p = [];
    fprintf('Skipped! [%s]\n', str);
    return;
end

% disp(r);

p.fname=r.fname;
p.argin.full=strtrim(r.argin);
r.stars=r.stars(find(r.stars~=' '));
p.argouttype.full=strtrim([ r.const ' ' r.basetype r.stars ]);
p.argouttype.const=r.const;
p.argouttype.basetype=r.basetype;
p.argouttype.stars=r.stars;

% parse input arguments
p.argin.args=regexp(p.argin.full,'(?<full>\w[^,]+[\w\*])','names');
p.argin.full='';
for i=1:numel(p.argin.args),

    r=regexp(p.argin.args(i).full,sprintf('(?<const>const)? ?(?<basetype>%s) ?(?<stars>[\\* ]*) ?(?<argname>\\w*)',types),'names');
    if isempty(r) || length(r) < 1
        p = [];
        fprintf('Skipped-II! [%s]\n', str);
        return;
    end
    
    r.stars=r.stars(find(r.stars~=' '));
    p.argin.args(i).type.full=strtrim([ r.const ' ' r.basetype r.stars ]);
    p.argin.args(i).full=strtrim([ r.const ' ' r.basetype r.stars ' ' r.argname ]);
    p.argin.args(i).type.const=r.const;
    p.argin.args(i).type.basetype=r.basetype;
    p.argin.args(i).type.stars=r.stars;
    p.argin.args(i).argname=r.argname;

    % reassemble full declaration
    p.argin.full=[ p.argin.full p.argin.args(i).full ];
    if i<numel(p.argin.args),
        p.argin.full=[ p.argin.full ', ' ];
    end
    
end

% reassemble full declaration
p.full=strtrim([ p.argouttype.full ' ' p.fname '(' p.argin.full ')' ]);

return
