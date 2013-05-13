function tcpip_calc_server(port)
% tcpip_calc_server  - Waits for "feval" calls from remote machine.
% Listens on port 26789 and services tcpip_feval calls from the
% remote macine. This server script uses tcpip_feval_server(...)
%  
  if nargin<1,
    port=26789;
  end
  sock=tcpip_servsocket(port);
  if sock==-1,
    disp(sprintf('Cant use local socket no:%d as server',port));
    return;
  end
  try,
    while 1,
      serv=-1;
      while serv==-1,
	serv=tcpip_listen(sock);
	pause(0.1);
      end
      
      if serv==-1,
	disp(sprintf('Cant get connection on socket on port no:%d as server',port));
	tcpip_close all;
	return;
      end
      
      while 1,
      
	disp 'Waiting for next remote call...'
  
	str='';
	while length(str)<20,
	  str=tcpip_viewbuff(serv,20);
	  if tcpip_status(serv)==0,	break; end
	  pause(0.1);
	end
	if tcpip_status(serv)==0,	break; end
	disp '... somethink is coming... recive and execute...'
	while 1,
	  tcpip_feval_server(serv);
	end
      end
      tcpip_close(serv);
      disp 'Closed current connection!!'
    end
  catch
    disp 'Error...? Close ALL tcpip connections'
    tcpip_close all;
  end
  
    
