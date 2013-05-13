function newStr=ReplaceLineTerminators(str, newTerminator)
% newStr=ReplaceLineTerminators(str, newTerminator)
%
% Replace either mac, windows, or linux style line terminators in the
% string "str" with the the specified line terminator and return the
% result. Strings containing Macintosh, Windows or Unix line terminators
% (see below), or any mix thereof are handled as input.
%
% The argument "newTerminator" may be either the terminator characters
% which you want to substitute in or else a string identifying a 
% platform. If you provide terminators, these must be a combination of [
% \f\n\r\t\v], while terminators for the desired platform can be specified
% by the identifiers in teh middle column below:
%
%       Macintosh:      'Mac', 'Macintosh', 'OS 9' 'OS9'        CR      13
%       Windows:        'Win', 'Windows', 'DOS', 'MSDOS'        CRLF    13 10
%       Unix:           'Unix','Linux', 'BSD', 'OS X', 'OSX'    LF      10
%
% SEE ALSO: BreakLines

% HISTORY
%
% 12/09/03  awi     Wrote it.
% 10/06/05  awi     Note here cosmetic changes by dgp made between 12/9/03 and 10/6/05.
% 22/10/11  dcn     Updated help and now supporting arbitrary terminators
%

mac=1;
win=2;
linux=3;        % "unix" is a MATLAB function so we use "linux" instead.  

platforms(mac).index=mac;
platforms(mac).name='Macintosh';
platforms(mac).break=char(13);
platforms(mac).aliases=upper({platforms(1).name, 'mac', 'os9', 'os 9', 'cr', platforms(1).break});

platforms(win).index=win;
platforms(win).name='Windows';
platforms(win).break=char([13 10]);
platforms(win).aliases=upper({platforms(2).name, 'win', 'dos', 'msdos', 'crlf', platforms(2).break});

platforms(linux).index=linux;
platforms(linux).name='Linux';
platforms(linux).break=char(10);
platforms(linux).aliases=upper({platforms(3).name, 'unix', 'bsd', 'lf', platforms(3).break});

% find the desired terminator from the newTerminator argument.
platformIndex=0;
for i=mac:linux
    if any(strcmpi(newTerminator,platforms(i).aliases))
        platformIndex=i;
        desiredTerminator = platforms(i).break;
        break;
    end
end
if platformIndex==0
    % input can be the desired terminators, which must be some combination
    % of whitespace characters [ \f\n\r\t\v]. Check that it is.
    nonWhiteSpace = regexp(newTerminator,'\S','once');
    if ~isempty(nonWhiteSpace)
        error('Unrecognized "newTerminator" argument value');
    else
        desiredTerminator = newTerminator;
    end
end

% if we started out as either Mac, Window, or Linux the result will be
% either Mac or Linux 
strMacOrLinux=strrep(str,platforms(win).break, platforms(mac).break);
% if we are given either Mac or Linux the result will be Linux.
strLinux=strrep(strMacOrLinux, platforms(mac).break, platforms(linux).break);

%now all the eol characters in the string are linux eols so replace
%them with the requested values.
newStr=strrep(strLinux, platforms(linux).break, desiredTerminator);
