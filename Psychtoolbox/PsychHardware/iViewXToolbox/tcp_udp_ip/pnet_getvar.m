function varargout=pnet_getvar(fid)
% PNET_GETVAR - Gets any matlab variable (number-, cell-, struct- or object -array) from pnet.
%
% Syntax:  
%  
%   [variable1, variable2,.....]=pnet_getvar(con)
%  
%   Receives matlab variables over a TCP connection that was sent with
%   PNET_PUTVAR. This variable transfer uses its own non standard protocol.
%   Can also be used with UDP packets if the variable/variables fits inside
%   an UDP packet.
%  
% Se also:  PNET_PUTVAR, PNET, PNET_REMOTE
%  


%
%   This file(s) is part of the tcp_udp_ip toolbox (C) Peter Rydesäter et al.  
%   et al.  1998-2003 for running in MATLAB(R) as scripts and/or plug-ins.
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%
%   In addition, as a SPECIAL EXCEPTION, Peter Rydesäter, SWEDEN,
%   gives permission to link the code of this program with any library,
%   and distribute linked combinations of it. You must obey the GNU
%   General Public License in all respects for all of the code in the
%   tcp_udp_ip toolbox and any . If you modify any source file included,
%   you may extend this exception to your version of the file, but you are
%   not obligated to do so.  If you do not wish to do so, delete this exception
%   statement from your version. This exception makes it possible to use
%   pnet.c (.dll) as a plug-in as it is intended and let it be (dynamical)
%   linked to MATLAB(R) or any compiled stand alone application.
%  
  
  if nargout~=1,
    for n=1:nargout,
      varargout{n}=pnet_getvar(fid);
    end
    return;
  end
  while 1,
    dataclass=pnet(fid,'readline',1024);
    switch dataclass,
     case {'double' 'char' 'int8' 'int16' 'int32' 'uint8' 'uint16' 'uint32'}
      datadims=double(pnet(fid,'Read',1,'uint32'));
      datasize=double(pnet(fid,'Read',datadims,'uint32'));
      VAR=pnet(fid,'Read',datasize,dataclass);
      break;
     case '--matfile--'
      tmpfile=[tempname,'.mat'];
      VAR=[];
      try,
	bytes=double(pnet(fid,'Read',[1 1],'uint32'));
	pnet(fid,'ReadToFile',tmpfile,bytes);
	load(tmpfile);
      end
      try,
	delete(tmpfile);
      end
      break;
     otherwise
    end
  end
  varargout{1}=VAR;
  return;
  
