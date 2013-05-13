function popmail_demo(site,user,pass)
% popmail_demo - Demo that read mail from pop mail server (not delete).
%                The first lines of each mail will be printed out.
%
% Syntax:
%    popmail_demo(site,user,password)
%  or
%    popmail_demo site user password
%  or
%    popmail_demo
%
%  In the last case you will be asked for intput parameters.
%
%  Version: 2002-02-01 Uppgraded to use the API in tcpiptoolbox 2.x, from tcpiptoolbox 1.x API
%
  
  if nargin<3,
    site=input('Input adress to pop server:','s');
    user=input('user:','s');
    pass=input('password: (sorry, will be displayed on the screen)','s');
  end
  % CONNECT
  fid=pnet('tcpconnect',site,110);
  if fid==-1,
    disp 'Cant connect to server!';
    return;
  end
  %LOGIN
  read_mresp(fid);
  pnet(fid,'printf','USER %s\n',user);
  read_mresp(fid);
  pnet(fid,'printf','PASS %s\n',pass);
  read_mresp(fid);
  pnet(fid,'printf','STAT\n');
  pnet(fid,'readline');
  all=0;
  % READ HEADERS OF FIRST 50 mail.
  for a=1:50 ,
    if all, break; end
    pnet(fid,'printf','TOP %d 0\n',a);
    s='';
    b=0;
    while strncmp(s,'.',1)==0,
      s=pnet(fid,'readline');
      b=b+1;
      if strncmp(s,'-ERR',4),
        all=1;
        nummes=a-1;
        break;
      elseif strncmp('Subject:',s,8);
        subject{a}=s(9:end);
      elseif strncmp('From:',s,5);
        from{a}=s(6:end);
      end
    end
    mlines(a)=b-1;
  end
  disp(sprintf('Number of mails in mailbox on server: %d\n',nummes));
  % DISPLAY HEADER AND READ MESSAGE LINES
  for a=1:nummes,
    disp '###############################################################'
    disp(sprintf('Subject: %s\n   From: %s',subject{a},from{a}));
    disp '---------------------------------------------------------------'
    %Request next line in current mail message.
    pnet(fid,'printf','RETR %d\n',a);
    for b=1:mlines(a),
      s=pnet(fid,'readline');
    end
    b=0;
    atflag=0;   % Set "not displaying" flag to false.
    while strcmp(s,'.')==0,
      b=b+1;
      s=pnet(fid,'readline');
      
      % Don't display attachment lines and similar stuff.
      if atflag==0 && b<100,
        disp(s);
      end
      if strncmp(s,'Content-',8),
        atflag=1;  % If start of attachment set "not displaying" flag.
      end
    end
  end
  pnet(fid,'close');
  return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read and check that response is OK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok=read_mresp(fid)
  s=pnet(fid,'readline');
  if strncmp(s,'+OK',3)==0,
    pnet(fid,'close');
    error('Response error!');
  end
  return;
