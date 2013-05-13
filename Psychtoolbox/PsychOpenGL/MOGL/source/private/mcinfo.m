function [ M, C ] = mcinfo( funcp )

% MCINFO  From a parsed C declaration, extract information needed to
%         write M-file wrapper and C interface function
% 
% [ M, C ] = mcinfo( funcp )

% 24-Jan-2005 -- created;  adapted from code in autocode.m (RFM)
% 06-Feb-2007 -- Modified; can now handle OpenAL as well. (MK)
% 24-Mar-2011 -- Modified; Perform 64-bit safe wrapping of encoded memory pointers. (MK)
% 28-Aug-2012 -- Modified; Handle datatypes GLint64, GLuint64 and GLsync. (MK)

% make lowercase function name, e.g., glVertex3dv --> gl_vertex3dv
k=min(find(funcp.fname~=lower(funcp.fname))); %#ok<MXFND>
lowerfn=lower([ funcp.fname(1:k-1) '_' funcp.fname(k:end) ]);

% initialize properties of C function
C.interfacename=lowerfn;    % name of interface function
C.arg_in={};                % input arguments
C.arg_out={};               % return argument
C.hasptrreturn = 0;         % 1 = Return argument is encoded memory pointer.

% initialize proerties of M-file wrapper
M.fname=funcp.fname;        % M-file filename
M.arg_in={};                % input arguments
M.arg_in_check={};          % code to check input arguments
M.arg_out={};               % return arguments
M.arg_out_init={};          % code to initialize return arguments
M.mogl_in={};               % input arguments in call to moglcore
M.mogl_out={};              % return arguments in call to moglcore
M.allocate=0;		        % flag indicating whether arguments to
							% moglcore.mexmac must be allocated manually

% if return type of function declaration is not void, then add return
% arguments to C function and M-file
if strcmp(funcp.argouttype.basetype,'void')==0,
	% scalar return argument
	if isempty(funcp.argouttype.stars) && strcmp(funcp.argouttype.basetype,'GLsync')==0,
		C.arg_out{1}=sprintf('plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);\n\t*mxGetPr(plhs[0])=(double)');
		M.arg_out{1}='r';
		M.mogl_out{1}='r';
	% pointer return argument
    else
        C.hasptrreturn = 1;
		% Old style: 32-bit pointer encoded in 32 bit uint32: C.arg_out{1}=sprintf('plhs[0]=mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);\n\t*(unsigned int *)mxGetData(plhs[0])=(unsigned int)');
        
        % New style: 32 or 64 bit pointer encoded in double (which can hold
        % up to 64-bit if properly mangled):
        % TODO: We could also return a mxUINT64_CLASS type aka uint64(), which
        % would be nicer, but only after we decide to completely drop
        % support for Matlab versions < R2007a.
		C.arg_out{1}=sprintf('plhs[0]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);\n\t*(double *)mxGetData(plhs[0])=PsychPtrToDouble((void*) ');
		M.arg_out{1}='r';
		M.mogl_out{1}='r';
	end
else
	C.arg_out{1}='';
end

% step through input arguments
for j=1:numel(funcp.argin.args),

	% skip void input arguments
	if strcmp(funcp.argin.args(j).type.full,'void'),
		continue
	end

	% provide an argument name, if none specified in declaration
	if isempty(funcp.argin.args(j).argname),
		funcp.argin.args(j).argname=sprintf('arg%d',j);
	end

    % 64-Bit (unsigned) integer scalar input argument?
	if isempty(funcp.argin.args(j).type.stars) && (strcmp(funcp.argin.args(j).type.basetype,'GLint64') || strcmp(funcp.argin.args(j).type.basetype,'GLuint64')),
        % These need special treatment, as a double scalar can't fully
        % represent a 64-Bit integer without loss of precision:
		C.arg_in{end+1}=sprintf('(%s) *((%s*) mxGetData(prhs[%d]))',funcp.argin.args(j).type.full,funcp.argin.args(j).type.basetype,j-1);
		M.arg_in{end+1}=funcp.argin.args(j).argname;
        if strcmp(funcp.argin.args(j).type.basetype,'GLint64')
            % Wrapper needs to cast to int64()
            M.mogl_in{end+1} = sprintf('int64(%s)', funcp.argin.args(j).argname);
        else
            % Wrapper needs to cast to uint64()
            M.mogl_in{end+1} = sprintf('uint64(%s)', funcp.argin.args(j).argname);
        end
        
    % GLsync handle scalar input argument?
	elseif isempty(funcp.argin.args(j).type.stars) && strcmp(funcp.argin.args(j).type.basetype,'GLsync'),
        % New style: Memory pointers are encoded opaque inside double
        % scalar values. This is 64-bit safe.
        % TODO: We could also use a mxUINT64_CLASS type aka uint64(), which
        % would be nicer, but only after we decide to completely drop
        % support for Matlab versions < R2007a.
        C.arg_in{end+1}=sprintf('(%s) PsychDoubleToPtr(mxGetScalar(prhs[%d]))',funcp.argin.args(j).type.full,j-1);
		M.arg_in{end+1}=funcp.argin.args(j).argname;
		M.arg_in_check{end+1}=sprintf('if ~strcmp(class(%s),''double''),\n\terror([ ''argument ''''%s'''' must be a pointer coded as type double '' ]);\nend\n',funcp.argin.args(j).argname,funcp.argin.args(j).argname);
		M.mogl_in{end+1}=funcp.argin.args(j).argname;
        
	% other scalar input argument?
    elseif isempty(funcp.argin.args(j).type.stars),
		C.arg_in{end+1}=sprintf('(%s)mxGetScalar(prhs[%d])',funcp.argin.args(j).type.full,j-1);
		M.arg_in{end+1}=funcp.argin.args(j).argname;
		M.mogl_in{end+1}=funcp.argin.args(j).argname;

	% honourary scalar?  (i.e., pointer to a GLU struct)
	elseif ~isempty(strfind(funcp.argin.args(j).type.basetype,'GLU')),
		% Old style: Only for 32 bit systems: C.arg_in{end+1}=sprintf('(%s)(unsigned int)mxGetScalar(prhs[%d])',funcp.argin.args(j).type.full,j-1);
        % New style: Memory pointers are encoded opaque inside double
        % scalar values. This is 64-bit safe.
        % TODO: We could also use a mxUINT64_CLASS type aka uint64(), which
        % would be nicer, but only after we decide to completely drop
        % support for Matlab versions < R2007a.
        C.arg_in{end+1}=sprintf('(%s) PsychDoubleToPtr(mxGetScalar(prhs[%d]))',funcp.argin.args(j).type.full,j-1);
		M.arg_in{end+1}=funcp.argin.args(j).argname;
		M.arg_in_check{end+1}=sprintf('if ~strcmp(class(%s),''double''),\n\terror([ ''argument ''''%s'''' must be a pointer coded as type double '' ]);\nend\n',funcp.argin.args(j).argname,funcp.argin.args(j).argname);
		M.mogl_in{end+1}=funcp.argin.args(j).argname;

	% pointer?
	else

		% decide how to cast input argument, if necessary
		switch funcp.argin.args(j).type.basetype,
			case { 'GLdouble' 'GLclampd' 'ALdouble' 'ALclampd' }
				mcast='double';
			case { 'GLfloat' 'GLclampf' 'ALfloat' 'ALclampf' }
				mcast='single';
			case { 'GLint64' }
                mcast='int64';
			case { 'GLuint64' }
                mcast='uint64';
            case { 'GLint' 'GLsizei' 'ALint' 'ALsizei' }
				mcast='int32';
			case { 'GLuint' 'GLenum' 'GLbitfield' 'GLUnurbs' 'GLUtesselator' 'GLUquadric' 'ALuint' 'ALenum' 'ALbitfield' }
				mcast='uint32';
			case { 'GLshort' 'ALshort' }
				mcast='int16';
			case { 'GLushort' 'ALushort'}
				mcast='uint16';
			case { 'GLbyte' 'ALbyte' }
				mcast='int8';
			case { 'GLubyte' 'GLchar' 'GLboolean' 'ALubyte' 'ALchar' 'ALboolean' }
				mcast='uint8';
			case { 'GLvoid' 'void' 'ALvoid' }
				mcast='';
            case { 'GLhandleARB' 'GLhandle' 'ALhandleARB' 'ALhandle' 'GLsync'}
                % We use a double as an opaque data type for these, as they
                % can be pointers or integers of 32-bit or 64-bit size.
                % Safe approach is to use double, as it can contain 64-bits
                % in principle. However it is not clear how sane it is to
                % pass these data types to the scripting environment at
                % all...
                mcast='double';
                fprintf('Warning: double cast for GLhandle style type. Doublecheck results for %s\n', funcp.fname);
			otherwise
				error(sprintf('don''t know how to handle pointer type %s',funcp.argin.args(j).type.full)); %#ok<SPERR>
		end

		% add input argument to C function
		C.arg_in{end+1}=sprintf('(%s)mxGetData(prhs[%d])',funcp.argin.args(j).type.full,j-1);

		% check wither pointer is const;  if not, assume data will be
		% written to this input argument, so initialize it
		if isempty(funcp.argin.args(j).type.const),
			thisallocate=1;
			M.allocate=1;
			M.arg_out{end+1}=funcp.argin.args(j).argname;
			M.arg_out_init{end+1}=sprintf('%s = %s(0);',funcp.argin.args(j).argname,mcast);
		else
			thisallocate=0;
			M.arg_in{end+1}=funcp.argin.args(j).argname;
		end

		% if pointer is not void, and is const, then typecast argument
		if ~isempty(mcast) && ~thisallocate,
			M.mogl_in{end+1}=sprintf('%s(%s)',mcast,funcp.argin.args(j).argname);
		else
			M.mogl_in{end+1}=sprintf('%s',funcp.argin.args(j).argname);
		end

	end

end

return
