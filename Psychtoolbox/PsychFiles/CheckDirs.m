function CheckDirs(dirs,mode)
% CheckDirs(dirs,mode):
% Iterates over all fields in struct 'dirs' and checks whether the Contents are existing directory addresses.
% mode == 'check': Display an error if this is not the case.
% mode == 'make' : Creates the directory if it doesnt exist.
% Default mode is 'check'.
% JvR 2008-03-31

psychassert(isstruct(dirs),'Piss off with your type ''%s''',class(dirs)) % Input must be a struct.

if nargin==1
    mode = 'check';                                            	% No mode input? Use default (check)
end

switch mode
    case 'make'
        structfun(@makedirs,dirs);                              % Apply function makedirs to each field in the struct
    case 'check'
        structfun(@checkem,dirs);                               % Apply function checkem to each field in the struct
    otherwise
        error('\nUnknown MODE input argument "%s"! \n\nRecognised inputs: ''check'' or ''make''.',mode);
end
end

function checkem(in)
psychassert(isdir(in),'Directory %s not found',in);          % Directory doesn't exist; error.
end

function makedirs(in)
if ~isdir(in)
    mkdir(in);                                              % Directory doesn't exist; create.
    disp(sprintf('Directory %s created.',in));
end
end
