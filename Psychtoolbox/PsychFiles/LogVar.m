function fname = LogVar(var,name,dirnm)
% fname = LogVar(var,name,dirnm)
%
% Turns a variable VAR into a string that evaluates back to the original
% variable, e.g. when feeding it to eval()
% NAME is the name of the variable that will be used in this string
% String will be saved in a text file in the directory DIRNM with as
% filename the name of the variable NAME and the time at which the file is
% written.
%
% DN 2008

str         = Var2Str(var,name);

logtijd     = clock;
fname       = [name ' log ' StrPad(logtijd(1),4,0) '-' StrPad(logtijd(2),2,0) '-' StrPad(logtijd(3),2,0) ' ' StrPad(logtijd(4),2,0) StrPad(logtijd(5),2,0) '.txt'];

fid         = fopen(fullfile(dirnm,fname), 'wt');
fprintf(fid,'%s',str);
fclose(fid);
