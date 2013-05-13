function tcpip_feval_server(ip_fid)
% tcpip_feval_server(ip_fid)  - Waits for "feval" calls on connection.
% Listens on ip_fid for remote call packets send from the client
% with   tcpip_feval(....)
%  

  callvar=tcpip_getvar(ip_fid);
  switch callvar{1},
   case 0,
    feval(callvar{2:end});
   case 1,
    [a1]=feval(callvar{2:end});
    ret={a1};
   case 2,
    [a1,a2]=feval(callvar{2:end});
    ret={a1,a2};
   case 3,
    [a1,a2,a3]=feval(callvar{2:end});
    ret={a1,a2,a3};
   case 4,
    [a1,a2,a3,a4]=feval(callvar{2:end});
    ret={a1,a2,a3,a4};
   case 5,
    [a1,a2,a3,a4,a5]=feval(callvar{2:end});
    ret={a1,a2,a3,a4,a5};
   case 6,
    [a1,a2,a3,a4,a5,a6]=feval(callvar{2:end});
    ret={a1,a2,a3,a4,a5,a6};
    tcpip_sendvar(ip_fid,ret);
   otherwise,
    [a1,a2,a3,a4,a5,a6]=feval(callvar{2:end});
    ret={a1,a2,a3,a4,a5,a6};
  end
  if callvar{1}>0,
    tcpip_sendvar(ip_fid,ret);
  end
  
  
