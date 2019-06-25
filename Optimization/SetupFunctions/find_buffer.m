function gen = find_buffer(gen,subnet,horizon)
%% calculate the upper and lower buffer thresholds for each storage system
n_g = length(gen);
gen_names = cell(n_g,1);
network_names = fieldnames(subnet);
for i = 1:1:n_g
    gen_names(i,1) = {gen(i).Name};
end
for net = 1:1:length(network_names)
    %identify equipment at each subNet node
    for m = 1:1:length(subnet.(network_names{net}).nodes)
        equip = subnet.(network_names{net}).Equipment{m};
        for j = 1:1:length(equip)
            if isfield(gen(equip(j)).QPform,'Stor') && any(strcmp(gen(equip(j)).QPform.states,'U'))%storage with a buffer state
                if isfield(gen(equip(j)).VariableStruct,'Buffer')
                    buff_perc = gen(equip(j)).VariableStruct.Buffer;% percentage for buffer on storage
                else
                    gen(equip(j)).VariableStruct.Buffer = 0;
                    buff_perc = 0;
                end
                if strcmp(network_names{net},'Hydro')
                    %hydro
                    dischargeCapacity = (gen(equip(j)).VariableStruct.MaxGenFlow + gen(equip(j)).VariableStruct.MaxSpillFlow)/12.1; %flow rate in 1000 ft^3 converted to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                    dischargeCapacity = dischargeCapacity*horizon/10; %amount the resevoir can discharge in 10% of the dispatch horizon.
                    buffer = min((buff_perc/100)*gen(equip(j)).QPform.Stor.UsableSize,dischargeCapacity);
                else
                    buffer = (buff_perc/100)*gen(equip(j)).QPform.Stor.UsableSize;
                end
                gen(equip(j)).QPform.link.bineq(end-1) = -buffer; %lower buffer ineq :  -SOC - W <= -Buffer becomes W>= buffer -SOC
                gen(equip(j)).QPform.link.bineq(end) = gen(equip(j)).QPform.Stor.UsableSize-buffer; %upper buffer ineq :  SOC - Z <= (UB-Buffer)
                gen(equip(j)).QPform.U.ub = buffer;
                gen(equip(j)).QPform.L.ub = buffer;
            elseif isfield(gen(equip(j)).QPform,'Stor')
                gen(equip(j)).VariableStruct.Buffer = 0;
            end
        end
    end
end