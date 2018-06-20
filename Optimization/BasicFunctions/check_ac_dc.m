function [gen,n_g] = check_ac_dc(gen,buildings,test_data)
%double check & add AC_DC conversion if neccessary
dir = strrep(which('check_ac_dc.m'),fullfile('Optimization','BasicFunctions','check_ac_dc.m'),'');
has_dc = false;
has_ac = false;
has_ac_dc = false;
need_ac = false;
need_dc = false;

n_g = length(gen);
for i = 1:1:n_g
    if strcmp(gen(i).Type,'AC_DC') 
        has_ac_dc = true;
    end
    if isfield(gen(i).Output,'DirectCurrent')
        has_dc = true;
    end
    if isfield(gen(i).Output,'Electricity') || strcmp(gen(i).Type,'Hydro Storage')
        has_ac = true;
    end
end
if isfield(test_data,'Demand') && isfield(test_data.Demand,'E')
    need_ac = true;    
end
if isfield(test_data,'Demand') && isfield(test_data.Demand,'DC')
    need_dc = true;    
end
if ~isempty(buildings)
    need_ac = true;
    n_b = length(buildings);
    for i = 1:1:n_b
        if isfield(buildings(i),'DCloads')
            need_dc = true;
        end
    end
end
if ~has_ac_dc && (need_ac && ~has_ac) || (need_dc && ~has_dc)
    load(fullfile(dir,'System Library','AC_DC','IdealConverter.mat'))
    gen(n_g+1) = component;
    n_g = n_g+1;
end
end%ends function check_ac_dc