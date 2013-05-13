function varargout = Explode(vect,splitvect,mode)
% outputcell = Explode(vect,splitvect,mode)
% [out1,out2,...] = Explode(vect,splitvect,mode)
% 
% Split a string or vector VECT by a string, scalar or vector SPLITSVECT.
% Output OUTPUTCELL will be returned as an array of cells.
% 
% Both VECT and SPLITSVECT must be a string or numeric vector of length 1
% or more.
% 
% The way this function handles its input is different for strings and
% numeric datatypes.
%
% strings
% SPLITSVECT can be specified using two formats (MODE):
% If MODE is not specified, 'fixed' will be used
%
%   'fixed'    : a fixed string, which will be used literally to split
%                VECT. Some character combinations have a special meaning:
%                   \b Backspace
%                   \f Form feed
%                   \n New line
%                   \r Carriage return
%                   \t Horizontal tab
%                example:
%                   s1 = 'resp32|me too|"get over here"';
%                   s2 = Explode(s1,'|'); or s2 = Explode(s1,'|','fixed');
%                   s2(:)
%                   ans = 
%                       'resp32'
%                       'me too'
%                       '"get over here"'
%                
%   'variable' : if you want to split the input with on a string that
%                matches a certain pattern.
%                specify this pattern using a regular expression,
%                see "doc regexp"
%                example:
%                   s1 = '|fff|ja|fdf|er|fft|fofr|';
%                   s2 = Explode(s1,'f.f','variable'); % . means any character
%                   s2(:)
%                   ans = 
%                       '|'
%                       '|ja|'
%                       '|er|fft|'
%                       'r|'
%
% numerical vectors
% SPLITSVECT can be specified using two formats (MODE):
% If MODE is not specified, 'fixed' will be used
% NB: splitting floats works, but precision will be lost
%
%   'fixed'    : a fixed number or sequence of numbers, which will be used
%                literally to split the vector.
%                example:
%                   n1 = [inf 33 12 45 13 nan 46 74 12 45 15 64];
%                   n2 = Explode(n1,[12 45]); of n2 = Explode(s1,[12 45],'fixed');
%                   n2(:)
%                   ans = 
%                       [Inf    33]
%                       [ 13   NaN    46    74]
%                       [ 15    64]
%                   n2 = Explode(n1,12);
%                   n2(:)
%                   ans = 
%                       [Inf    33]
%                       [ 45    13   NaN    46    74]
%                       [ 45    15    64]
% 
%   'variable' : If you want to split on a sequence of numbers that
%                contains (a) wildcard(s). Each wilcard can match a single
%                item in the vector.
%                NaN specifies a wildcard. When using 'variable' for
%                numeric input, you can thus no longer use NaN in your
%                pattern.
%                example:
%                   n1 = [inf 33 12 45 13 nan 46 74 12 35 64 13 21];
%                   n2 = Explode(n1,[12 nan 13],'variable');
%                   n2(:)
%                   ans = 
%                       [Inf    33]
%                       [NaN    46    74    12    35    64    13    21]
%
%                   n2 = Explode(n1,[12 nan nan 13],'variable');
%                   n2(:)
%                   ans = 
%                       [Inf    33    12    45    13   NaN    46    74]
%                       [21]

% DN 2008
% DN 2008-07-31 Added checking of MATLAB version
% DN 2010-12-29 Fixed version check now that subversion is over .9...
% DN 2011-09-05 Fixed a lot of numeric splitting. Now deals with floats
%               correctly (though precision is still lost), supports
%               variable number of outputs and clarified meaning of
%               wildcard for number splitting.

if ~IsOctave
    mver = ver('matlab');
    psychassert(str2double(mver.Version(1))>=7 && str2double(mver.Version(3:end))>2,'Need at least MATLAB R2006b');
else
    % we don't have the 'split' option in regexp, at least in version 3.2.4
    % strsplit only splits on single delimiters, so not an option
    error('Explode is not supported on Octave');
end

if nargin==3
    switch mode
        case 'fixed'
            qfixed   = true;
        case 'variable'
            qfixed   = false;
        otherwise
            error('Invalid mode "%s" specified, use "fixed" (default) or "variable"',mode);
    end
else
    qfixed   = true;
end

if ischar(vect)
    if qfixed
        % prepare splitsvect for use with regexp - e.g. convert return to
        % the regexp pattern for return
        splitvect = regexprep(splitvect,'\\([^bfnrt])','\\\\$1');      % escape backslash unless followed by a special character
    end
    strc = regexp(vect,splitvect,'split');
    
elseif isnumeric(vect)
    dtype       = class(vect);
    vect        = num2str(vect(:)');                                    % make sure its a columnvector and convert to string
    qisNaN      = all(isnan(splitvect(:)));
    splitvect   = num2str(splitvect(:)');
    splitvect   = regexprep(splitvect,'\s+','\\s+');
    if ~qfixed
        psychassert(~qisNaN,'When mode==''variable'', split vector cannot contain only nan, cannot split on only a wildcard')
        splitvect   = regexprep(splitvect,'[nN]a[nN]','[\\w\\d\\.]+?');
    end
    
    split       = regexp(vect,['(?<=\s)' splitvect '(?=\s)'],'split');
    func        = @(x)cast(str2num(x),dtype); %#ok<ST2NM>
    strc        = cellfun(func,split,'UniformOutput',false);
else
    error('input type %s not supported',class(vect));
end

if nargout <= 1
    varargout{1} = strc;
else
    varargout = strc;
end
