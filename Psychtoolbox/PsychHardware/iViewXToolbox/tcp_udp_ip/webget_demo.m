function [pagedata,hdr]=webget_demo(url)
% WEBGET_DEMO - Demo webpage downloader, returns a webpage as a string
%
% Syntax:
%    webget_demo url_adress
%  or
%    str=webget_demo('url_adress')
%  
%  Version: 2002-02-01 for the tcpiptoolbox 2.x API
%
if nargin==0,    url='http://www.mathworks.com'; end
MAXPAGESIZE=1024*1024;
%Decode URL
if strncmp(url,'http:',5), url=url(6:end); end
if strncmp(url,'//',2),    url=url(3:end); end
[host,page]=strtok(url,'/');
[host,port]=strtok(host,':');
if length(port)>1, port=port(2:end); else port='80'; end
if length(page)==0, page='/'; end
%Connect to web server
con=pnet('tcpconnect',host,port);
if con==-1, error 'Bad url or server down.....'; end
disp(['Connected to: ' host]);
if 1,
  pnet(con,'setwritetimeout',1);
  pnet(con,'setreadtimeout',30);
  %Send request
  pnet(con,'printf','GET %s HTTP1.0\n\n',page);
  pnet(con,'read',MAXPAGESIZE,'view');   % Read data to local buffer with 30 sec timout.
  pnet(con,'setreadtimeout',0.1); %Change timeout
  hdrstr=''; hdr='';
  %Remove header
  while length(hdrstr) || length(hdr)==0,
    stat=pnet(con,'status');
    hdrstr=pnet(con,'readline');
    hdr=sprintf('%s%s\n',hdr,hdrstr);
  end
  %Return page
  pnet(con,'setreadtimeout',60);
  pagedata=pnet(con,'read',MAXPAGESIZE);
end
pnet(con,'close');
