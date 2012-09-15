function [retcode] = sendmail(from, to, subj, mesg, varargin)
% SENDMAIL Send Internet e-mail
%   Using SENDMAIL (which needs the TCP/UDP/IP toolbox, freely available from
%   http://petrydpc.ite.mh.se/tools/) it is possible to send e-mail messages
%   from the Matlab prompt. This can be useful when you want to be notified
%   when large jobs terminate on remote machines.
%
% Usage
%   [RETCODE] = SENDMAIL(FROM, TO, SUBJ, MESG [, FILENAME]) send an e-mail with
%   subject SUBJ and contents MESG to the email address in TO, with
%   the From: field set to FROM.
%
%   FROM, TO and SUBJ are Matlab strings. MESG can be either a
%   string, or a cell array of strings.
%
%   RETCODE is 0 if the message was sent succesfully, otherwise it
%   has the value -1.
%
%   The optional argument FILENAME is a string containing the
%   filename of the file to be attached to the message. It may
%   contain path information.
%
% Examples
%   sendmail('me@some.where.com','me@else.where.org','Job Finished!','It took 4 days, 12 hours,6 minutes.')
%
%   sendmail('me@some.where.com','me@else.where.org','Job Finished!','Finaly...', 'results.mat')
%
% Installation
%   1. Make sure the TCP/UDP/IP toolbox is somewhere in your PATH
%   2. Adjust the SMTPSERVER, SMTPSERVERPORT and CLIENTIP variables to reflect your situation.

% Author: Joris Portegies Zwart
%         http://www.science.uva.nl/~portegie/
% Version 3.1
% Date:   December 24, 2002

% Contributors: Simon Bridger

% Revisions :
% 3.1 Removed bug in filename extraction
% 3.0 Added attachment capabilities
% 2.1 SJB
%     - Add "Content-Transfer-Encoding" to work w/outlook97/exchange server
%     - Made server info optional
%     - Added possibility of mesg in cell format
% 2.0 Update to use the new TCP/UDP/IP Toolbox 2.0.2
% 1.3 Added server response checking
% 1.2 Changed <CR> to <CRLF>, as in RFC821
% 1.1 Corrected minor syntax errors
% 1.0 First version

% Copyright (C) 2002 Joris Portegies Zwart
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details:
%
%         http://www.gnu.org/copyleft/gpl.html
  

ncmesg = nargchk(4, 5, nargin);
if ~isempty(ncmesg)
  error(ncmesg)
end

SMTPSERVER = 'mail.server.com';       % IP adress of SMTP server
SMTPSERVERPORT = 25;                      % Port number, ususally 25
CLIENTIP = 'my.computerip.com';     % IP adress of client machine, usually not important
CRLF = [char(13) char(10)];

S = rand('state');
rand('state',sum(100*clock));
MIMEBREAK = char(floor(97 + 25 * rand(1, 25)));
rand('state', S);

if iscellstr(mesg) % Message is a cellarray of strings
  msg_string='';
  for i=1:length(mesg)
    msg_string=[msg_string, mesg{i}, CRLF];
  end
  mesg=msg_string;
end

if (nargin==5)
  fullname = varargin{1};
  fsid = findstr(fullname, filesep);
  if ~isempty(fsid)
    filename = fullname( (fsid(end)+1):end);
  else
    filename = fullname;
  end
  disp('Encoding attachment.')
  b64data = base64encode(fullname);
  hasattachment = 1;
else
  hasattachment = 0;
end

try
   smtp_connection = pnet('tcpconnect', SMTPSERVER, SMTPSERVERPORT);
   str = '';
   while (isempty(str))
      str=pnet(smtp_connection,'readline');
   end
   disp(['Connection with ' SMTPSERVER ':' num2str(SMTPSERVERPORT) ' established.']);
catch
   disp('Network Connection Failed, Aborting sendmail')
   retcode = -1;
   return
end


Header =...
    ['From: ' from CRLF ...
     'To: ' to CRLF ...
     'Subject: ' subj CRLF ...
     'MIME-Version: 1.0' CRLF ...
     'Content-type: multipart/mixed; boundary="' MIMEBREAK '"' CRLF];


Message =...
    ['This is a multi-part message in MIME format.' CRLF CRLF ...
     '--' MIMEBREAK CRLF ...
     'Content-type: text/plain; charset=us-ascii' CRLF ...
     'Content-transfer-encoding: 7bit' CRLF CRLF ...
     mesg CRLF CRLF];

if hasattachment
  Attachm =...
      ['--' MIMEBREAK CRLF ...
       'Content-type: application/octet-stream; name="' filename '"' CRLF ...
       'Content-transfer-encoding: base64' CRLF CRLF];
  for n = 1:length(b64data)
    Attachm = [Attachm b64data{n} CRLF];
  end
else
  Attachm = [];
end
   
Closing =[CRLF '--' MIMEBREAK '--' CRLF CRLF '.'];

try
   sendcommand(smtp_connection, ['HELO ' CLIENTIP]);
   sendcommand(smtp_connection, ['MAIL FROM:' from ]);
   sendcommand(smtp_connection, ['RCPT TO:' to ]);
   sendcommand(smtp_connection, ['DATA']);
   sendcommand(smtp_connection, [Header Message Attachm Closing]);
   sendcommand(smtp_connection, ['QUIT']);
catch
   disp('Aborting sendmail')
   pnet(smtp_connection,'close');
   retcode = -1;
   return
end

pnet(smtp_connection,'close');
disp('Email sent.')
retcode = 0;

return



% SENDCOMMAND %
function [str] = sendcommand(smtp_connection, command)

CRLF = [char(13) char(10)];
str = '';

%    This function call should return the number of character written,
%    but doesn't. The proper syntax would be
%
% numc_written = pnet(smtp_connection,'printf', '%s', [command CRLF]);
%
%    but we have to use

pnet(smtp_connection,'printf', '%s', [command CRLF]);

while (isempty(str))
   str = pnet(smtp_connection, 'readline');
end

%    Uncomment the next line to echo all server responses
% disp(str)

resp_top = str2num(str(1));
switch resp_top
	case {2, 3} % 2xx means succesful completion, 3xx means awaiting input
   	return
	case{4} % 4xx means transient failure, could try again?
   	disp(['Transient failure issuing command ' command])
   	disp(['The server said ' str])
   	error('ERROR!')
	case{5} % 5xx means permanent failure
   	disp(['Permanent failure issuing command ' command])
   	disp(['The server said ' str])
   	error('ERROR!')
	otherwise % Unspecified server behaviour
   	disp('Unknown server response, trying to continue ...');
   	return
end

return



% BASE64ENCODE %
function [out] =  base64encode(filename)
% This function reads a file as a bitstream, 76 groups of 6 bits at
% a time. Each group gets encoded as an 8-bit character from the
% base64 dictionary. Each line of 76 characters is stored in a cell
% of the cell array out.

BASE64 = ['ABCDEFGHIJKLMNOPQRSTUVWXYZ' ...  % The Base64 dicitionary.
	  'abcdefghijklmnopqrstuvwxyz' ...
	  '0123456789+/'];
BIT2DEC = 2.^[5 4 3 2 1 0];
NBITS = 6 * 76;                  % Number of bits per line

fid = fopen(filename, 'r', 'b'); % Open file as big-endian
status = fseek(fid, 0, 'eof');   % Set pointer to end of file
fsize = ftell(fid);              %  to determine size.
frewind(fid)                     % Reset pointer

numlines = ceil((4 * fsize/3)/76);
out = cell(1,numlines);

for n = 1:(numlines-1) % Encoding the first numlines-1 lines is easy.
  [fc, count] = fread(fid, NBITS, 'ubit1');
  out{n} = BASE64(1 + BIT2DEC * reshape(fc, 6, 76));
end

% The last line needs a little more work. First, read remaining bits:
[fc, count] = fread(fid, NBITS, 'ubit1');
%Pad last part of stream so it's length is multiple of 6 bits:
fc = [fc ; zeros(6*ceil(count/6) - count, 1)];
% Encode final group of 6 bits
out{numlines} = BASE64(1 + BIT2DEC * reshape(fc, 6, length(fc)/6));
% Pad encoded stream such with '=' characters number of characters
% in encoding is integral multiple of 4
numpad = rem(fsize * 8, 24)/8;    
out{numlines} = [out{end} repmat('=', 1, numpad)];

fclose(fid);
return
