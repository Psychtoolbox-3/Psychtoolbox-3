function fid=tcpip_servopen(port)
%TCPIP_SERVOPEN opens a blocking TCPIP server socket.
%
% Call with:
%
%   fid = tcpip_servopen(port)
%
% where:
%
%   port   is the port number of the socket you wish to open.
%
%   fid      is the retuned handle to the socket, if successful,
%            -1 if unsuccessful.
%
% This function is obsolete.  It has been replaced with TCPIP_SERVSOCKET
% and TCPIP_LISTEN.  It is only included for backward compatibility.
%
% TCPIP_SERVOPEN opens a TCPIP server socket with the requested PORT
% number.  This function blocks (waits) until a client socket has
% established a connection before returning FID.
%
% EXAMPLE:
%
%    disp('Waiting for someone to connect to port 4444...');
%    fid = tcpip_servopen(4444);
%    if fid = -1
%      error('TCP/IP Connection error!');
%    end
%    disp('Connection! :-)');
%    disp('Now sending a message and closing!')
%    tcpip_write('Hello!!',string(10)); %string(10) generates newline
%    tcpip_close(fid);
%
% See also:  TCPIP_SERVSPCKET, TCPIP_LISTEN, TCPIP_OPEN
%            TCPIP_CLOSE
%
  fid=-1;
  serv_fid = tcpip_servsocket(port);
  if serv_fid == -1
    return;
  end
  try,
    while fid == -1
      fid = tcpip_listen(serv_fid);
      pause(1);
    end
  catch,
    tcpip_close(serv_fid);
    return;
  end
  tcpip_close(serv_fid);
  return;
  
