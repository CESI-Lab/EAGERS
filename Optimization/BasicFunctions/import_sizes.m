function size_info = import_sizes(filename,z_names,l_names,multiplier)
[cat,heading,z_info] = read_eio(filename);

zone = nonzeros((1:length(cat))'.*strcmp('Zone Sizing Information',cat));
z_objects = z_info{zone};
sys = nonzeros((1:length(cat))'.*strcmp('System Sizing Information',cat));
s_objects = z_info{sys};
comp = nonzeros((1:length(cat))'.*strcmp('Component Sizing Information',cat));
c_objects = z_info{comp};

size_info.zone_flow = zeros(length(z_names),1);
size_info.zone.calc_load = zeros(length(z_names),2);
size_info.zone.user_load = zeros(length(z_names),2);
size_info.zone.calc_flow = zeros(length(z_names),2);
size_info.zone.user_flow = zeros(length(z_names),2);
k = nonzeros((1:length(heading{zone}))'.*strcmp('User Des Air Flow Rate {m3/s}',heading{zone}));
for i = 1:1:length(z_objects(:,1))
    z = nonzeros((1:length(z_names))'.*strcmpi(fix_name(z_objects{i,1}),z_names));
    new_flow = str2double(z_objects{i,k});
    size_info.zone_flow(z) = max(size_info.zone_flow(z),new_flow/multiplier(z));
    if strcmpi(z_objects{i,2},'heating')
        a = 1;
    else
        a = 2;
    end
    size_info.zone.calc_load(z,a) = str2double(z_objects{i,3});
    size_info.zone.user_load(z,a) = str2double(z_objects{i,4});
    size_info.zone.calc_flow(z,a) = str2double(z_objects{i,5});
    size_info.zone.user_flow(z,a) = str2double(z_objects{i,6});
    if isempty(z_objects{i,8})
        size_info.zone.hour_of_day(z,a) = 12;
    else
        size_info.zone.hour_of_day(z,a) = time_convert(z_objects{i,8});
    end
end

size_info.loop_flow = zeros(length(l_names),1);
size_info.loop.user_load = zeros(length(l_names),2);
size_info.loop.calc_flow = zeros(length(l_names),2);
size_info.loop.user_flow = zeros(length(l_names),2);
k = nonzeros((1:length(heading{sys}))'.*strcmp('User Des Air Flow Rate [m3/s]',heading{sys}));
for i = 1:1:length(s_objects(:,1))
    l = nonzeros((1:length(l_names))'.*strcmpi(s_objects{i,1},l_names));
    new_flow = str2double(s_objects{i,k});
    size_info.loop_flow(l) = max(size_info.loop_flow(l),new_flow);
    if strcmpi(s_objects{i,2},'heating')
        a = 1;
    elseif strcmpi(s_objects{i,2},'cooling')
        a = 2;
    end
    size_info.loop.user_load(l,a) = str2double(s_objects{i,4});
    size_info.loop.calc_flow(l,a) = str2double(s_objects{i,5});
    size_info.loop.user_flow(l,a) = str2double(s_objects{i,6});
    size_info.loop.hour_of_day(l,a) = time_convert(s_objects{i,8});
end

for i = 1:1:length(c_objects(:,1))
    size_info.component.type(i,1) = c_objects(i,1);
    size_info.component.name(i,1) = c_objects(i,2);
    size_info.component.field(i,1) = c_objects(i,3);
    size_info.component.value(i,1) = str2double(c_objects{i,4});
end
end%Ends import_sizes

function [cat,heading,z_info] = read_eio(filename)
if ~contains(filename,'.eio')
    filename = strcat(filename,'.eio');
end
text = fileread(filename);
n_line = strfind(text,char(10));
n_l = length(n_line);
lines = cell(n_l,1);
s = 1;
for i = 1:1:length(n_line)
    lines(i) = {text(s:n_line(i)-2)};
    s = n_line(i)+1;
end

j = 0;
for i = 2:1:length(lines)-1
    com = strfind(lines{i},',');
    if strcmp(lines{i}(1),'!')
        j = j+1;
        cat(j,1) = {lines{i}(4:com(1)-2)};
        if length(com)==1
            headings = {rmv_spaces(lines{i}(com(1)+1:end))};
        else
            headings = {rmv_spaces(lines{i}(com(1)+1:com(2)-1))};
            for ii = 2:1:length(com)-1
                headings = [headings;rmv_spaces(lines{i}(com(ii)+1:com(ii+1)-1))];
            end
            headings = [headings;rmv_spaces(lines{i}(com(end)+1:end))];
        end
        heading(j,1) = {headings};
        z_info(j,1) = {cell(0,length(com))};
    else
        cat_i = nonzeros((1:length(cat))'.*strcmp(rmv_spaces(lines{i}(1:com(1)-1)),cat));
        if isempty(cat_i)
            cat_i = length(cat);
        end
        k = length(z_info{cat_i}(:,1))+1;
        if length(com)==1
            z_info{cat_i,1}{k,1} = rmv_spaces(lines{i}(com(1)+1:end));
        else
            for ii = 1:1:length(com)-1
                z_info{cat_i,1}{k,ii} = rmv_spaces(lines{i}(com(ii)+1:com(ii+1)-1));
            end
            z_info{cat_i,1}{k,length(com)} = rmv_spaces(lines{i}(com(end)+1:end));
        end
    end
end
end%Ends function read_eio

function name = rmv_spaces(name)
while ~isempty(name) && strcmp(name(1),' ')
    name = name(2:end);
end
while ~isempty(name) && strcmp(name(end),' ')
    name = name(1:end-1);
end
end%Ends function rmv_spaces

function hod = time_convert(str)
col = strfind(str,':');
h = str2double(str(col(1)-3:col(1)-1));
m = str2double(str(col(1)+1:col(1)+2));
s = str2double(str(col(2)+1:end));
hod = h+m/60+s/3600;
end%Ends function time_convert

function name = fix_name(name)
name = rmv_spaces(name);
if ~isempty(name)
    name = strrep(strrep(name,' ','_'),',','');
    name = strrep(strrep(name,'-','_'),'/','_');
    name = strrep(strrep(name,':','_'),'#','number');
    name = strrep(strrep(name,'(',''),')','');
    if any(strcmp(name(1),{'1';'2';'3';'4';'5';'6';'7';'8';'9';'0';}))
        name = strcat('a_',name);
    end
    if length(name)>60
        name = name(1:60);
    end
end
end%Ends function fix_name