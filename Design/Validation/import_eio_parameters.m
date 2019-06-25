function param = import_eio_parameters(filename)
text = fileread(filename);
n_line = strfind(text,char(10));
lines = cell(length(n_line),1);
s = 1;
for i = 1:1:length(n_line)-3
    lines(i) = {text(s:n_line(i)-2)};
    s = n_line(i)+1;
    if i == 1
        %skip first line
    elseif strcmp(lines{i}(1),'!')
        brack1 = strfind(lines{i},'<');
        brack2 = strfind(lines{i},'>');
        title = fix_name(lines{i}(brack1+1:brack2-1));
        param.(title) = [];
        com = strfind(lines{i},',');
        nH = length(com);
        headings = cell(nH,1);
        for j = 1:1:nH-1
            headings(j) = {lines{i}(com(j)+1:com(j+1)-1)};
        end
        headings(nH) = {lines{i}(com(nH)+1:end)};
        for j = 1:1:nH
            name = headings{j};
            units = [];
            brack3 = strfind(name,'{');
            if ~isempty(brack3)
                brack4 = strfind(name,'}');
            else
            	brack3 = strfind(name,'(');
                if ~isempty(brack3)
                    brack4 = strfind(name,')');
                else
                    brack3 = strfind(name,'[');
                    if ~isempty(brack3)
                        brack4 = strfind(name,']');
                    end
                end
            end
            if ~isempty(brack3)
                units = name(brack3+1:brack4-1);
                if brack4<length(name)
                    name = [name(1:brack3-1),name(brack4+1:end)];
                else
                    name = name(1:brack3-1);
                end
            end
            name = fix_name(name);
            param.(title).(name) = [];
            if ~isempty(units)
                param.(title).(name).units = units;
            end
        end
    else
        com = strfind(lines{i},',');
        title = fix_name(lines{i}(1:com(1)-1));
        f_names = fieldnames(param.(title));
        for j=1:1:length(f_names)
            if isfield(param.(title).(f_names{j}),'value')
                k = length(param.(title).(f_names{j}).value)+1;
            else
                k = 1;
            end
            if j<length(com)
                val = lines{i}(com(j)+1:com(j+1)-1);
            else
                val = lines{i}(com(j)+1:end);
            end
            num = str2double(val);
            if isnan(num)
                if strcmp(val,'N/A')
                	param.(title).(f_names{j}).value(k,1) = num;
                else
                    param.(title).(f_names{j}).value(k,1) = {rmv_spaces(val)};
                end
            else
                param.(title).(f_names{j}).value(k,1) = num;
            end
        end            
    end
end

end%Ends function import_idf_parameters

function name = fix_name(name)
name = rmv_spaces(name);
name = strrep(strrep(name,'-','_'),'/','_');
name = strrep(strrep(strrep(name,' ','_'),':','_'),'#','number');
name = strrep(strrep(name,'(',''),')','');
if any(strcmp(name(1),{'1';'2';'3';'4';'5';'6';'7';'8';'9';'0';}))
    name = strcat('a_',name);
end
end%Ends function fix_name

function name = rmv_spaces(name)
while ~isempty(name) && strcmp(name(1),' ')
    name = name(2:end);
end
while ~isempty(name) && strcmp(name(end),' ')
    name = name(1:end-1);
end
end%Ends function rmv_spaces