function psychassert(varargin)
% psychassert(expression, ...) - Replacement for Matlab 7 builtin assert().
% This is hopefully useful for older Matlab installations and
% for the Octave port:
%
% If the assert-function is supported as builtin function on
% your Matlab installation, this function will call the builtin
% "real" assert function. Read the "help assert" for usage info.
%
% If your Matlab lacks a assert-function, this function
% will try to emulate the real assert function. The only known limitation
% of our own assert wrt. Matlabs assert is that it can't handle the MSG_ID
% parameter as 2nd argument. Passing only message strings or message
% formatting strings + variable number of arguments should work.
%

% History:
% 01/06/09 mk Wrote it. Based on the specification of assert in Matlab 7.3.

if exist('assert', 'builtin')==5
  % Call builtin implementation:
  builtin('assert', varargin{:});
else
  % Use our fallback-implementation:
  if nargin < 1
      error('Not enough input arguments.');
  else
      expression = varargin{1};
      if ~isscalar(expression) || ~islogical(expression)
          error('The condition input argument must be a scalar logical.');
      end
      
      % Expression true?
      if ~expression
          % Assertion failed:
          if nargin < 2
              error('Assertion failed.');
          else
              if nargin < 3
                  emsg = sprintf('%s\n', varargin{2});
                  error(emsg); %#ok<SPERR>
              else
                  emsg = sprintf(varargin{2}, varargin{3:end});
                  error(emsg); %#ok<SPERR>                  
              end
          end
      end
  end
end

return;
