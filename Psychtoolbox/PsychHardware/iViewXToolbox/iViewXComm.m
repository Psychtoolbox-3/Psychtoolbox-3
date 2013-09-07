function [result, ivx]=iViewXComm(cstr, ivx, sendstr)

% wrapper function for in and outgoing udp connections


% pnet messages
% /* Different status of a con_info struct handles a file descriptor    */
% #define STATUS_FREE       -1
% #define STATUS_NOCONNECT   0    // Disconnected pipe that is note closed
% #define STATUS_TCP_SOCKET  1
% #define STATUS_IO_OK       5    // Used for IS_... test
% #define STATUS_UDP_CLIENT  6
% #define STATUS_UDP_SERVER  8
% #define STATUS_CONNECT     10   // Used for IS_... test
% #define STATUS_TCP_CLIENT  11
% #define STATUS_TCP_SERVER  12
% #define STATUS_UDP_CLIENT_CONNECT 18
% #define STATUS_UDP_SERVER_CONNECT 19
%
%

result=-1;
switch ivx.connection
    case 1, % ethernet
        switch lower(cstr)
            case 'open',
                if isempty(ivx.udp)
                    % open (and close) udp connection for receiving data from iView computer
                    % we should split opening, checking, receiving and closing of course

                    % Open  udpsocket and bind udp port adress to it.
                    ivx.udp=pnet('udpsocket',ivx.localport);
                    pnet(ivx.udp,'setreadtimeout',ivx.udpreadtimeout);
                    if ~isempty(ivx.udp)
                        stat=pnet(ivx.udp,'status');
                        if stat>0
                            result=ivx.isconnected;
                        end
                    end
                end
            case 'receive',
                if ~isempty(ivx.udp)
                    % should we instead check status?
%                     stat=pnet(ivx.udp,'status');

                    % Wait/Read udp packet to read buffer
                    len=pnet(ivx.udp,'readpacket');
                    %len=pnet(udp,'readpacket',[],'noblock');
                    % if len>0 fprintf('Len: %d\n', len); end

                    if len>0,
%                         [ip,port]=pnet(ivx.udp,'gethost');
                        % if packet larger then 1 byte then read maximum of 1000 doubles in network byte order
                        %  data=pnet(udp,'read',ivx.udpmaxread,'double');
                        data=pnet(ivx.udp,'read',ivx.udpmaxread);
                        % data
                        result=data;
                    end
                end
            case 'receivelast',
                if ~isempty(ivx.udp)
                    % should we instead check status?
                    % stat=pnet(ivx.udp,'status');
                    
                    while 1
                        % Wait/Read udp packet to read buffer
                        len=pnet(ivx.udp,'readpacket');
                        
                        %len=pnet(udp,'readpacket',[],'noblock');
                        % if len>0 fprintf('Len: %d\n', len); end
                        
                        if len > 0
                            % [ip,port]=pnet(ivx.udp,'gethost');
                            % if packet larger then 1 byte then read
                            % maximum of 1000 doubles in network byte order
                            % data=pnet(udp,'read',ivx.udpmaxread,'double');
                            data=pnet(ivx.udp,'read',ivx.udpmaxread);
                            % data
                            result=data;
                        else
                            break
                        end
                    end
                end
            case 'close',
                if ~isempty(ivx.udp)

                    pnet(ivx.udp,'close');
                    ivx.udp=[];
                end
                result=1;
            case 'check',
                if ~isempty(ivx.udp)
                    stat=pnet(ivx.udp,'status');
                    if stat>0
                        result=ivx.isconnected;
                    end
                end
            case 'send',

                % open udp connection, and send command string
                % if this is too time consuming, we should open a port and keep it open
                % for the time of the experiment
                % code stolen from udp_send_demo.m
                udp=pnet('udpsocket',ivx.socket);
                if udp~=-1,
                    try, % Failsafe
                        pnet(udp,'udpconnect',ivx.host,ivx.port);
                        %         [ip,port]=pnet(udp,'gethost')
                        %         stat=pnet(udp,'status')
                        pnet(udp,'write',[sendstr char(10)]);        % Write to write buffer
                        pnet(udp,'writepacket',ivx.host,ivx.port);   % Send buffer as UDP packet

                    catch,
                        pnet('closeall');
                        disp(lasterr)
                    end
                    %     [ip,port]=pnet(udp,'gethost')
                    %     stat=pnet(udp,'status')
                    pnet(udp,'close');
                    result=1;
                end


            otherwise,
                error([mfilename ': unknown command for ethernet communication']);
        end
    case 2,
        % assume we use serial communication
        error([mfilename ': serial communication not implemented yet']);


    otherwise,
        error([mfilename ': unknown connection type']);
end

