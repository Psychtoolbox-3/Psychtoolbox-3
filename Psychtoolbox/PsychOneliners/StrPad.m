function str = StrPad(in,len,char)
% str = StrPad(in,length,padchar)
% (pre)pads IN with CHAR to sepcified length LEN. If inputs IN or PADCHAR
% are numerical, they will be converted to to string. If input is too long,
% it is truncated from the start to specified length.
%
% DN 2007

if isnumeric(in) && length(in)==1 && in==round(in)
    % convert to string
    in = num2str(in);
end
if isnumeric(char) && length(char)==1
    % convert to string
    char = num2str(char);
end
if ischar(in)
    % check that we have a string
    inlen = length(in);
    if inlen > len
        % truncate
        b = inlen - len; % string is b characters too long
        str = in(b+1:end);
    elseif inlen == len
        % string is right length already
        str = in;
    else
        % pre-pad
        b = len - inlen; % string is b characters too short
        str = [repmat(char,1,b) in];
    end
else
    error('input must be char or scalar integer');
end
    
