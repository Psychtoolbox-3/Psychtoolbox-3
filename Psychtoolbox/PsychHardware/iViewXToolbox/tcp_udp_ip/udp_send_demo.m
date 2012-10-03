function udp_send_demo(fun,host,port)
% UDP_SEND_DEMO - a demo that sends a squence of doubles in network byte order to local or remote host
%
% Syntax:
%   UDP_SEND_DEMO
% or
%   UDP_SEND_DEMO function_string
% or
%   UDP_SEND_DEMO function_string hostname
% or
%   UDP_SEND_DEMO function_string hostname portnumber
%
% Default values:
%    function_string  is by default sin(0:0.1:6) that will be evaluated and transmitted 
%                     as a sequence of network byte ordered doubles (or generated datatype)
%
%    hostname         is by default localhost but can be any hostname if you whant to send
%                     the packet to an other host.
%
%    portnumber       is by default 3333.
%
%
% The purpose of this demo is to illustrate how a udp packat can be created, filled with numbers
% and then transmitted to any host and udp port. Use this demo together with udp_plotter_demo
% that receives and plott the packets of numbers.
%  
% Example:
%
% udp_send_demo sin(0:0.1:50)./(0:0.1:50) plotterhost 33333
%

if nargin<1, fun='sin(0:0.1:6)'; end
if nargin<2, host='localhost'; end
if nargin<3, port='3333'; end

data=evalin('caller',fun);
udp=pnet('udpsocket',1111)
if udp~=-1,
  try, % Failsafe
    pnet(udp,'write',data);              % Write to write buffer
    pnet(udp,'writepacket',host,port);   % Send buffer as UDP packet
    udp
  end
  pnet(udp,'close');
end
