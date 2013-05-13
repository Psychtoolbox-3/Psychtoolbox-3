function www_sql()
%
% WWW-server as interface to ALIS data in MySQL
%
  tcpip_close all;
  while 1,
    sock=-1;
    while sock<0,
      disp 'Open socket 8080...'
      sock=tcpip_servsocket(8080);
      if(con<0), pause(0.1); end
    end
    while tcpip_status(sock),
      disp 'Listen for connections...'    
      con=-1;
      while con<0,
	con=tcpip_listen(sock);
	if(con<0), pause(0.1); end
      end
      disp ':-)   Connection....'
      try, www_read_do(con); end   %% Read from connection...
      tcpip_close(con);
    end
    tcpip_close all;
  end
  return;
  
%
% Read http data from connection and decide what to do.
%
function www_read_do(con)
  str='';
  start_time=now;
  while length(str)==0 && tcpip_status(con)>=10 && now-start_time < 10,
    str=tcpip_readln(con,1000);
  end
  if length(str)==0,
    disp 'Read error...'
    return;
  end
  
  if strncmp(upper(str),'GET /',5),
    www_write_main(con);
  end      
  return;

%
% Write main page!
%
function www_read_do(con)
  pg=sprintf(['HTTP/1.1 200 OK\n'...
	      'Content-Type: text/html\n\n'...
	      '<html>\n'...
	      '<h1>WEBSERVER DEMO</h1>'...
	      '</html>']);
  tcpip_write(con,pg);
  return;
  

