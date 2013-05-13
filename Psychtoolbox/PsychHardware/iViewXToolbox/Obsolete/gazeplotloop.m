host='10.0.1.2';
host=[]; % default
port=[]; %default 

while KbCheck; end

tEnd=GetSecs+20;

x=400;
y=400;

mag=8;

while ~KbCheck && Getsecs<tEnd
    x=x+round(-1*mag+rand*2*mag);
    y=y+round(-1*mag+rand*2*mag);
    
    if x<100, x=100; end
    if x>800, x=800; end

    if y<100, y=100; end
    if y>600, y=600; end

   
    datastr=['[' num2str(x) ' ' num2str(y) ']'];
    udp_send_demo_OSX(datastr,host,port);
    
%     WaitSecs(0.1);
end
