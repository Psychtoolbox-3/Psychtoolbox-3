function output = PackColorImage(red,green,blue)	
% output = PackColorImage(red,green,blue)
% 
% Take three image planes (same dimension m by n)
% and pack them into an m by n by 3 image.
% 
% Particularly useful for fixing old code
% that used SCREEN('PutColorImage',....)
%
% 11/24/02  jmh, dhb  Wrote it.
% 06/12/12        dn  Use cat() for this

output = cat(3,red,green,blue);
