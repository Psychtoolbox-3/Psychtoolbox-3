function vec = Struct2Vect(struc,fieldnm)
% vec = Struct2Vect(struc,fieldnm)
% Traverses array of structs STRUC and returns data from all fields FIELDNM
% in vector VEC.
% Returns array vector if field contains numeric scalars or scalar structs,
% cell vector otherwise. If concatenating struct arrays, each struct must
% contain the same fields
%
% Example: for field data.field:
%   data(1).field = 23;
%   data(2).field = 56;
%   vec = Struct2Vect(data,'field')
%   vec = [23 56];
% 
% Example: for field data2.test:
%   data2(1).test = 23;
%   data2(2).test = 'd';
%   vec = Struct2Vect(data2,'test')
%   vec = {23 'd};
% 
% Example: for field data3.cells:
%   data3(1).cells = 23;
%   data3(2).cells = [42 45];
%   vec = Struct2Vect(data3,'cells')
%   vec = {23, [42 45]};
% 
% Example: for field data4.struc:
%   data4(1).struc.a = 24;
%   data4(2).struc.a = [42 45];
%   vec = Struct2Vect(data4,'struc')
%   vec(1).a = 24;
%   vec(2).a = [42 45];


% DN 2007
% DN 2008-07-30 Fixed handling of numeric vectors
% DN 2012-06-12 Better tests for numeric scalar, now tests all elements in
%               array. Can now handle struct arrays

isnumer = arrayfun(@(x) isnumeric(x.(fieldnm)),struc); isnumer = all(isnumer(:));
isstruc = arrayfun(@(x) isstruct(x.(fieldnm)),struc);  isstruc = all(isstruc(:));
isscala = arrayfun(@(x) isscalar(x.(fieldnm)),struc);  isscala = all(isscala(:));

if (isnumer || isstruc) && isscala
    vec = [struc.(fieldnm)];
else % not numeric scalar, make cellarray
    vec = {struc.(fieldnm)};
end
