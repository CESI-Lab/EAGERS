function costs = defaultCosts(costs,gen)
%DEFAULTCOSTS
if isempty(costs) || ~isfield(costs,'DiscountRate')
    costs.DiscountRate = 2;% implies a 2% discount rate for NPC calculations
end
n_g = length(gen);
if ~isfield(costs,'Equipment') || length(costs.Equipment)~=n_g
    costs.Equipment = [];
    for i = 1:1:n_g
        cost_per_kw = [];
        costs.Equipment(i).Name = gen(i).Name;
        switch gen(i).Type
            case {'CHP Generator';'Electric Generator';'Hydrogen Generator';}
                if isfield(gen(i).Output,'DirectCurrent')
                    cost_per_kw = 3000;
                    o_and_m = 100;
                else
                    cost_per_kw = 1000;
                    o_and_m = 50;
                end
            case 'Electrolyzer'
                cost_per_kw = 1000;
                o_and_m = 50;
            case 'Chiller'
                if strcmp(gen(i).Source,'Electricity')
                    % electric chiller
                    cost_per_kw = 100;
                    o_and_m = 10;
                else
                    % absorption chiller
                    cost_per_kw = 200;
                    o_and_m = 20;
                end
            case 'Heater'
                cost_per_kw = 100;
                o_and_m = 3;
            case 'Hydro Storage'
                cost_per_kw = 50;
                o_and_m = 5;
            case 'Hydrogen Storage'
                cost_per_kw = 50;
                o_and_m = 5;
            case 'Thermal Storage'
                cost_per_kw = 20;
                o_and_m = 1;
            case 'Electric Storage'
                cost_per_kw = 500;
                o_and_m = 10;
            case 'Solar'
                cost_per_kw = 500;
                o_and_m = 5;
        end
        if ~isempty(cost_per_kw)
            costs.Equipment(i).Cost = cost_per_kw*gen(i).Size;
            costs.Equipment(i).OandM = o_and_m*gen(i).Size;
            costs.Equipment(i).Financed = 100;
            costs.Equipment(i).LoanRate = 6;
            costs.Equipment(i).LoanTerm = 15;
        end
    end
end