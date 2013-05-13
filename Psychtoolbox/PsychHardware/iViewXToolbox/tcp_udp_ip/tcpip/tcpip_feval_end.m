function [varargout]=tcpip_feval_end(ip_fid)
% [arg1,arg2,...]=tcpip_feval_end(ip_fid) get return args from remote feval
% You must get the number of arguments you specified when calling
% tcpip_feval(...). If you specified 0 arguments dont call this function
% at all!
%
  ret=tcpip_getvar(ip_fid);
  varargout=ret;
  return;
  
  
