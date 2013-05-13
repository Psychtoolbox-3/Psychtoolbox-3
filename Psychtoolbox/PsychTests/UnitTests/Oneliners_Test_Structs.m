function success = Oneliners_Test_Structs

success = true;

% create some test variables
a.test = 4;
a.yar = 'test';
b.test = 'collision';
c.der  = 42;
d.test = 4;
d.c    = {3,'true',false};
e.c    = {3,'true',false};
notstruct = 'nope, not a struct';
f.r     = 4;
f(3).t  = 5;
f(4).e  = [];
g.r     = 4;
g.e     = 34;
g(2).r  = 5;
g(3).r  = 4;
data(1).field = 23;
data(2).field = 56;
data2(1).test = 23;
data2(2).test = 'd';
data3(1).cells = 23;
data3(2).cells = [42 45];
data4(1).struc.a = 24;
data4(2).struc.a = [42 45];

try
    % test: AddStructs
    qCaught = false;
    try
        % this should generate an error
        out = AddStructs(a,b);
    catch
        qCaught = true;
    end
    if ~qCaught
        success = false;
        fprintf('AddStructs didn''t flag collision when merging two structs that contain the same field\n');
    end
    clear out, clear qCaught
    
    qCaught = false;
    try
        % this should generate an error
        out = AddStructs(a,notstruct);
    catch
        qCaught = true;
    end
    if ~qCaught
        success = false;
        fprintf('AddStructs didn''t flag that one of the inputs was not a struct\n');
    end
    clear out, clear qCaught
    
    out = AddStructs(a,c);
    fields = fieldnames(out);
    if length(fields)>3 || ~any(strcmp('test',fields)) || ~any(strcmp('yar',fields)) || ~any(strcmp('der',fields))
        success = false;
        fprintf('AddStructs fieldnames in result not correct\ngot fields:\n');
        fprintf('  %s\n',fields{:});
    elseif out.test~=4 || ~strcmp(out.yar,'test') || out.der~=42
        success = false;
        fprintf('AddStructs values in fields not correct\n');
    end
    clear out, clear fields
    
    % test: AreStructsEqualOnFields
    out = AreStructsEqualOnFields(a,b,'test');
    if out
        success = false;
        fprintf('AreStructsEqualOnFields these two structs should not compare equal on the field ''test''\n');
    end
    out = AreStructsEqualOnFields(a,d,'test');
    if ~out
        success = false;
        fprintf('AreStructsEqualOnFields these two structs should compare equal on the field ''test''\n');
    end
    out = AreStructsEqualOnFields(d,e,'c');
    if ~out
        success = false;
        fprintf('AreStructsEqualOnFields these two structs should compare equal on the field ''test''\n');
    end
    clear out
    
    % test: CleanStruct
    out = CleanStruct(f);
    if ~(length(out)==2 && out(1).r==4 && isempty(out(1).t) && out(2).t==5 && isempty(out(2).r))
        success = false;
        fprintf('CleanStruct did not remove empty structs from struct array\n');
    elseif isfield(out,'e')
        success = false;
        fprintf('CleanStruct did not remove empty field from struct array\n');
    end
    clear out
    
    % test: FillEmptyFields
    out = FillEmptyFields(f,'-');
    if ~(all(strcmp({out.e},'-')) && out(2).r=='-' && out(3).r=='-')  % not going to test all, if these work rest will too
        success = false;
        fprintf('FillEmptyFields did not fill in empty fields in struct array correctly\n');
    end
    out = FillEmptyFields({'a',[],4,[]},'-');
    if ~(out{1}=='a' && out{2}=='-' && out{3}==4 && out{4}=='-')
        success = false;
        fprintf('FillEmptyFields did not fill in empty elements in cell array correctly\n');
    end
    clear out
    
    % test: GroupStructArrayByFields
    out = GroupStructArrayByFields(g,'r');
    if ~(length(out{1})==2 && isscalar(out{2}))
        success = false;
        fprintf('GroupStructArrayByFields did not group structs correctly\n');
    end
    clear out
    
    % test: Struct2Vect
    out = Struct2Vect(data,'field');
    if ~(isnumeric(out) && all(out==[23 56]))
        success = false;
        fprintf('Struct2Vect did not concatenate two scalar numerics correctly\n');
    end
    out = Struct2Vect(data2,'test');
    if ~(iscell(out) && out{1}==23 && out{2}=='d')
        success = false;
        fprintf('Struct2Vect did not concatenate a scalar numeric and a char correctly\n');
    end
    out = Struct2Vect(data3,'cells');
    if ~(iscell(out) && out{1}==23 && all(out{2}==[42 45]))
        success = false;
        fprintf('Struct2Vect did not concatenate a scalar numeric and numeric array correctly\n');
    end
    out = Struct2Vect(data4,'struc');
    if ~(isstruct(out) && out(1).a==24 && all(out(2).a==[42 45]))
        success = false;
        fprintf('Struct2Vect did not concatenate scalar structs correctly\n');
    end
    clear out
    
    
catch
    success = false;
    fprintf('error ocurred: "%s"\n',lasterr);
    
end
