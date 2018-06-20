function gen = max_utility_sellback(gen,subnet,test_data)
%identify upper bound for utility states (helps with scaling)
n_g = length(gen);
max_dc_gen = 0;
if isfield(subnet,'DirectCurrent')
    for n = 1:1:length(subnet.DirectCurrent.nodes)
        equip = subnet.DirectCurrent.Equipment{n};
        for j = 1:1:length(equip)
            if ~strcmp(gen(equip(j)).Type,'AC_DC')
                max_dc_gen = max_dc_gen+ max(0,gen(equip(j)).Size*gen(equip(j)).QPform.output.DC(1));
            end
        end
    end
end
max_ac_gen = 0;
if isfield(subnet,'Electrical')
    for n = 1:1:length(subnet.Electrical.nodes)
        equip = subnet.Electrical.Equipment{n};
        for j = 1:1:length(equip)
            if ~strcmp(gen(equip(j)).Type,'AC_DC')
                max_ac_gen = max_ac_gen+ max(0,gen(equip(j)).Size*gen(equip(j)).QPform.output.E(1));
            end
        end
    end
end
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Utility') && ~isempty(gen(i).QPform.states)%%avoid things like gas utility with no states
        %identify the network
        out = char(fieldnames(gen(i).QPform.output));
        if strcmp(out,'E')
            net = 'Electrical';
        elseif strcmp(out,'H')
            net = 'DistrictHeat';
        elseif strcmp(out,'C')
            net = 'DistrictCool';
        elseif strcmp(out,'W')
            net = 'Hydro';
        end
        if isfield(test_data,'Demand') && isfield(test_data.Demand,out)
            gen(i).QPform.X.ub = 10*max(sum(test_data.Demand.(out),2));%max Purchase
        else
            gen(i).QPform.X.ub = 1e6; %arbitrary upper bound that is not inf
        end
        if strcmp(net,'Electrical') && isfield(gen(i).QPform,'Y')
            gen(i).QPform.Y.ub = max_ac_gen + max_dc_gen;
        end
    end
    if strcmp(gen(i).Type,'AC_DC')
        gen(i).QPform.A.ub = max_ac_gen;
        gen(i).QPform.B.ub = max_dc_gen;
    end
end 