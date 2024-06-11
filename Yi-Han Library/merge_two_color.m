function output_struct = merge_two_color(struct1, name1, strunt2, name2)

if size(struct1,1)>size(strunt2,1)
    output_struct = struct(struct1);
    anoth_struct = struct(strunt2);
    newnameNumSML = strcat('NumSMLLabel', name1);
    newnameSML = strcat('SMLLabel', name1);
    anot_NumSML = strcat('NumSMLLabel', name2);
    anot_SML = strcat('SMLLabel', name2);
else
    output_struct = struct(struct2);
    anoth_struct = struct(strunt1);
    newnameNumSML = strcat('NumSMLLabel', name2);
    newnameSML = strcat('SMLLabel', name2);
    anot_NumSML = strcat('NumSMLLabel', name1);
    anot_SML = strcat('SMLLabel', name1);
end


[output_struct.(newnameNumSML)] = output_struct.NumSMLLabel1;
[output_struct.(newnameSML)] = output_struct.SMLLabel1;
for iiii = 1:size(output_struct,1)
    for jjjj = 1:size(anoth_struct,1)
        if anoth_struct(jjjj).partnum == output_struct(iiii).partnum
            output_struct(iiii).(anot_NumSML) = anoth_struct(jjjj).NumSMLLabel1;
            a(iiii).(anot_SML) = anoth_struct(jjjj).SMLLabel;
        end
    end
end
output_struct = rmfield(output_struct,'NumSMLLabel1');
output_struct = rmfield(output_struct,'SMLLabel');

