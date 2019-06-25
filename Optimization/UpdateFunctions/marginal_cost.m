function margin_cost = marginal_cost(gen,gen_disp,date)
%% This function estimates the marginal cost for each addittional kW of generating capacity
%%There are too many factors, e.g. storage, startup, co-production of heat, 
%%reserve capcity called upon at previous step..., to make a complete 
%%determination, so this is purely an estimate. The magin cost is sorted by
%%the cheapest available kW for generators that are already on, then
%%turning on generators one at a time in the cheapest order.
n = 4;%break into 4 segments
scale_cost = update_cost(date,gen); 
locked = gen_disp~=0;
n_s = length(date)-1;
margin_cost.Timestamp = date;
dt = (date(2:end) - date(1:end-1))*24;
s= {'E';'H';'C'};
s2 = {{'CHP Generator';'Electric Generator';};{'Heater'};{'Chiller';}};
n_g = length(gen);
dc = false;
for i = 1:1:n_g
    if (strcmp(gen(i).Type,'CHP Generator') || strcmp(gen(i).Type,'Electric Generator') || strcmp(gen(i).Type,'Hydrogen Generator')) && isfield(gen(i).VariableStruct.Startup,'DirectCurrent')
        dc = true;
    end
end
if dc
    s(end+1) = {'DC'};
    s2(end+1) = {{'CHP Generator';'Electric Generator';}};
end
max_out = gen_limit(gen,gen_disp,locked,dt);
for k = 1:1:length(s)
    inc = false(n_g,1);
    for i = 1:1:n_g
        if ismember(gen(i).Type,s2{k})
            inc(i) = true;
        end
    end
    if any(inc)
        out = s{k};
        margin_cost.(out).Capacity.SpinReserve = zeros(n_g,n_s,n);
        margin_cost.(out).Capacity.NonSpin = zeros(n_g,n_s,n);
        margin_cost.(out).Capacity.DemandResponse = zeros(n_g,n_s,n);
        margin_cost.(out).Cost.SpinReserve = zeros(n_g,n_s,n);
        margin_cost.(out).Cost.NonSpin = zeros(n_g,n_s,n);
        margin_cost.(out).Cost.DemandResponse = zeros(n_g,n_s,n);
        for i = 1:1:n_g
            if inc(i)
                s3 = fieldnames(gen(i).Output);
                output = s3{2}; %1st field is 'Capacity', 2nd field is primary output
                sr = (max_out.(s{k})(2:end,i)- gen_disp(2:end,i))';
                for t = 1:1:n_s
                    if sr(t)>1e-4
                        Range = linspace(gen_disp(t+1,i),max_out.(s{k})(t+1,i),n+1);
                        margin_cost.(out).Capacity.SpinReserve(i,t,:) = sr(t)/n;
                        cost = scale_cost(t+1,i)*Range./interp1(gen(i).Output.Capacity*gen(i).Size,gen(i).Output.(output),Range);
                        margin_cost.(out).Cost.SpinReserve(i,t,:) = (cost(2:n+1)-cost(1:n))/(sr(t)/n);%marginal cost per kW
                        margin_cost.(out).Cost.NonSpin(i,t,:) = inf;
                    else
                        Range = linspace(gen(i).VariableStruct.Startup.(output)(end),gen(i).Size,n);
                        cost = [0 scale_cost(t+1,i)*Range./interp1(gen(i).Output.Capacity*gen(i).Size,gen(i).Output.(output),Range)];
                        if isfield(gen(i).VariableStruct,'StartCost')
                            cost(2) = cost(2) + gen(i).VariableStruct.StartCost;
                        end
                        if isfield(gen(i).QPform,'constCost')
                            cost(2) = cost(2) + gen(i).QPform.constCost*scale_cost(t+1,i);
                        end
                        if strcmp(gen(i).Type,'Chiller') && isfield(gen(i).QPform.constDemand,'E')
                            cost(2) = cost(2) + gen(i).QPform.constDemand.E*min(nonzeros(margin_cost.E.Cost.SpinReserve(:,t,1)));
                        end
                        margin_cost.(out).Capacity.NonSpin(i,t,:) = Range - [0, Range(1:n-1)];
                        margin_cost.(out).Cost.NonSpin(i,t,:) = (cost(2:n+1)-cost(1:n))'./squeeze(margin_cost.(out).Capacity.NonSpin(i,t,:));%marginal cost per kW
                        margin_cost.(out).Cost.SpinReserve(i,t,:) = inf;
                    end
                end
            end
        end

    end 
end
%%make corrections for storage?
    
    %%find cost per kW of demand response
    
end %End function marginal_cost