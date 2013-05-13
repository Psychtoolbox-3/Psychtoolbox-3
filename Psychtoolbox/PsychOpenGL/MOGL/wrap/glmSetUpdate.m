function nproc = glmSetUpdate( newstate )

% glmSetUpdate  Set or check status of UNIX update process
% 
% usage:  nproc = glmSetUpdate( newstate )

% 09-Dec-2005 -- created (RFM)

% ---protected---

% get full path of setupdate.sh
cmd=which('setupdate.sh');

% if there's an input argument, start or stop the update process
if nargin>0,

    % check argument
    if ~isreal(newstate) || (newstate~=0 && newstate~=1),
        error('usage:  glmSetUpdate( [ 0 1 ] )');
    end

    % run script
    fullcmd=[ cmd ' ' num2str(newstate) ];
    if unix(fullcmd),
        warning('unable to execute ''%s''',fullcmd);
    end
    
end

% if there's a return argument, or no input arguments, then check number of
% update processes running
if nargout>0 || nargin==0,
    
    % check number of 'update' processes running
    [status,stdout]=unix(cmd);
    if status~=0,
        warning('unable to execute ''%s''',cmd);
        nproc=-1;
    end
    count=regexp(stdout,'^([0-9]*).*','tokens','once');
    nproc=str2num(count{1});

end

return
