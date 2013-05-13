function out = CleanStruct(temp,qrec)
% out = CleanStruct(temp,qrec)
%
% Deletes all empty structs from a struct array
% Deletes all empty fields from a struct (array)
% 
% if QREC == true (default is false), recurses through the struct and also
% cleans any nested structs it encounters in the struct
%
% DN    2007

if nargin==1
    qrec = false;
end

out = [];                                                           % initialize output struct 
c   = 1;                                                            % counter for number of non-empty elements in struct
if isstruct(temp)
    f       = fieldnames(temp);                                     % geeft een cellarray met daarin de veldnamen van de struct
    done    = false;                                                % boolean indicating whether a non-empty element has been found and processed

    for p=1:length(temp)                                            % nr of elements/levels in struct: struct(p).field
        for q=1:length(f)                                           % nr of fields in struct, for each field: struct(p).field{q}
            if ~isempty(temp(p).(f{q}))                             % non-empty element found, add it to output struct
                done = true;
                if isstruct(temp(p).(f{q})) && qrec                 % deeper level struct found, delete empty elements from that struct and add it to output struct
                    out(c).(f{q}) = CleanStruct(temp(p).(f{q}));
                else                                                % not a struct, add the element to output struct
                    out(c).(f{q}) = temp(p).(f{q});
                end
            end
        end
        if done
            % if a non-empty element has been found in this loop, this will be
            % true. Reset and increase counter by one for processing next
            % element.
            done = false;
            c = c + 1;
        end
    end
else
    out=temp;
end
