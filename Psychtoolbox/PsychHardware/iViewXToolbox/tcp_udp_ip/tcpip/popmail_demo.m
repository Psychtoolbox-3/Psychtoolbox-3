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
  
  nl=char(10);   % Define newline charecter variable.
  
  if nargin<3,
    site=input('Input adress to pop server:','s');
    user=input('user:','s');
    pass=input('password: (sorry will be shown)','s');
  end
  % CONNECT
  fid=tcpip_open(site,110);
  if fid==-1,
    disp 'Cant connect to server!';
    return;
  end
  %LOGIN
  if read_mresp(fid)==0, return; end
  tcpip_write(fid,['USER ',user,nl]);
  if read_mresp(fid)==0, return; end
  tcpip_write(fid,['PASS ',pass,nl]);
  if read_mresp(fid)==0, return; end
  tcpip_write(fid,['STAT',nl]);
  s=block_readln(fid);
  all=0;
  % READ HEADERS OF FIRST 50 mail.
  for a=1:50 ,
    if all, break; end
    tcpip_write(fid,['TOP ',num2str(a),' 0',nl]);
    s='';
    b=0;
    while strcmp(s,'.')==0,
      s=block_readln(fid);
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
  % DISPLAY HEADER AND READ MESSAGE LINES
  for a=1:nummes,
    disp '###############################################################'
    disp(sprintf('Subject: %s\n   From: %s',subject{a},from{a}));
    disp '---------------------------------------------------------------'
    %Request next line in current mail message.
    tcpip_write(fid,['RETR ',num2str(a),nl]);
    s='';
    
    for b=1:mlines(a),
      s=block_readln(fid);
    end
    b=0;
    atflag=0;   % Set "not displaying" flag to false.
    while strcmp(s,'.')==0,
      b=b+1;
      s=block_readln(fid);
      
      % Don't display attachment lines and similar stuff.
      if atflag==0 && b<100,
        disp(s);
      end
      if strncmp(s,'Content-',8),
        atflag=1;  % If start of attachment set "not displaying" flag.
      end
    end
  end
  tcpip_close(fid);
  return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read and check that response is OK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok=read_mresp(fid)
  ok=0;
  s=block_readln(fid);
  if strncmp(s,'+OK',3),
    ok=1;
  else
    disp('Response error!');
    tcpip_close(fid);
  end
  return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Busy wait loop function for reading of text line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=block_readln(fid)
  s='';
  while isempty(s),
    s=tcpip_readln(fid,10000);
  end
  % disp(['Response>>',s]);
  return;
