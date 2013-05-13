function tcpip_feval(ip_fid,varargin)
% tcpip_feval(ip_fid,no_ret_args,'fun_name',arg1,arg2,...) remote "feval"
% Sends to server request to execute 'fun_name' with argument....
% no_ret_args specifies number of return arguments. If return arguments is
% more then 0 then use tcpip_feval_end to get them later.
%
  tcpip_sendvar(ip_fid,varargin);
  
