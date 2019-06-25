function [verified,cost,binary_comb,disp_comb] = test_min_cases(qp,gen,options,net_demand,marginal,scale_cost,first_profile,dt)
%test cases with fewest # of generators then sort other by estimated cost
parallel = false;
if license('test','Distrib_Computing_Toolbox') 
    parallel = true;
end
n_g = length(gen);
stor_abbrev = cell(1,n_g);
stor_index = [];
for i = 1:1:n_g
    if isfield(gen(i).QPform,'Stor')
        f = fieldnames(gen(i).QPform.output);
        stor_abbrev(i) = f(1);
        stor_index(end+1) = i;
    end
end
c_red = best_eff(gen,scale_cost,marginal);    
[combinations,n_gen] = build_cases(qp,gen,options,net_demand,first_profile,dt);
c_remain = length(combinations(:,1));
verified = false(c_remain,1);
cost = nan(c_remain,1);
binary_comb = combinations>0;
disp_comb = zeros(c_remain,length(qp.organize));
if ~isempty(n_gen)
    test_more = true;
    best_cost = inf;
    l_t = 0;
    while test_more
        c_with_n_gen = nnz(n_gen(l_t+1:end)<=max(n_gen(1)+1,n_gen(l_t+1)));
        if parallel %&& c_with_n_gen>10
            need_test = [];
            for k = l_t+1:l_t+c_with_n_gen
                if isnan(cost(k)) || cost(k)<best_cost
                    need_test(end+1) = k;
                end
            end
            n_test = length(need_test);
            cost_par = cost(need_test);
            disp_comb_par = zeros(n_test,length(qp.organize));
            n_stor = length(stor_index);
            stor_pen_par = cell(n_test,1);
            c_test = combinations(need_test,:);
            parfor par_i = 1:n_test
                qp_test_par = disable_generators_step(qp,c_test(par_i,:)>0);%Disable generators here
                [x_par, flag1] = call_solver(qp_test_par);
                stor_pen_par{par_i} = zeros(2,n_g);
                if flag1 == 1
                    stor_pen_par{par_i} = zeros(2,n_g);
                    cost_par(par_i) = 0.5*x_par'*qp_test_par.H*x_par + x_par'*qp_test_par.f + sum(qp_test_par.constCost.*(c_test(par_i,:)>0));
                    disp_comb_par(par_i,:) = sort_solution_step(x_par,qp_test_par);
                    for j = 1:1:n_stor
                        s = qp_test_par.organize{stor_index(j)};
                        stor_pen_par{par_i}(1,stor_index(j)) = x_par(s);
                        stor_pen_par{par_i}(2,stor_index(j)) = qp_test_par.H(s,s);
                    end
                else
                    cost_par(par_i) = inf;
                end
            end
            cost(need_test) = cost_par;
            verified(need_test) = true;
            disp_comb(need_test,:) = disp_comb_par;
            best_cost = min(cost(1:l_t+c_with_n_gen));
            for nt = 1:1:n_test
                if ~isinf(cost_par(nt))
                    cost = cost_estimate(c_red,gen,cost,verified,marginal,stor_pen_par{nt},stor_abbrev,combinations,need_test(nt));%best case scenario for unverified cases
                end
            end
        else
            for k = l_t+1:1:l_t+c_with_n_gen
                if isnan(cost(k)) || cost(k)<best_cost
                    qp_test = disable_generators_step(qp,combinations(k,:)>0);%Disable generators here
                    [x, flag1] = call_solver(qp_test);
                    if flag1 == 1
                        cost(k) = 0.5*x'*qp_test.H*x + x'*qp_test.f + sum(qp_test.constCost.*(combinations(k,:)>0));
                        if cost(k)<best_cost
                            best_cost = cost(k);
                        end
                        stor_pen = zeros(2,n_g);
                        for j = 1:1:length(stor_index)
                            s = qp_test.organize{stor_index(j)};
                            stor_pen(1,stor_index(j)) = x(s);
                            stor_pen(2,stor_index(j)) = qp_test.H(s,s);
                        end
                        disp_comb(k,:) = sort_solution_step(x,qp_test);
                    else
                        cost(k) = inf;
                    end
                    verified(k) = true;
                    cost = cost_estimate(c_red,gen,cost,verified,marginal,stor_pen,stor_abbrev,combinations,k);%best case scenario for unverified cases
                end
            end
        end
        if l_t+c_with_n_gen == c_remain
            test_more = false;
        elseif ~any(isnan(cost)) && ~any(cost(l_t+c_with_n_gen+1:end)<best_cost)
            test_more = false;
        end
        l_t = l_t + c_with_n_gen;
    end
end
end%ends function test_minimum cases

function [combinations,n_gen] = build_cases(qp,gen,options,net_demand,first_profile,dt)
nn_list = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nn_abrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
include = {'Electric Generator';'CHP Generator';'Hydrogen Generator';'Heater';'Electrolyzer';'Chiller';'Absorption Chiller';'Cooling Tower';};%no utilities or storage systems are dispatchable
n_g = length(gen);
inc = zeros(n_g,1);
ninc = zeros(n_g,1);
for i = 1:1:n_g
    if ismember(gen(i).Type,include) && qp.Organize.Dispatchable(i) && qp.Organize.Enabled(i)
        inc(i) = i;
    elseif qp.Organize.Enabled(i)
        ninc(i) = i;
    end
end
inc = nonzeros(inc);
ninc = nonzeros(ninc);
n = length(inc);
combinations = zeros(2^n,n_g); %K is a matrix of all the possible generator on/off combinations 
if ~isempty(ninc)
    for j = 1:1:length(ninc)
        combinations(:,ninc(j)) = ninc(j); % all systems that are always included
    end
end

for j = 1:1:n %all combinations of generators are listed below
    z = 2^(n-j);
    r=0;
    while r+z<=2^n
        combinations(r+1:r+z,inc(j)) = inc(j);
        r=r+2*z;
    end
end
%sort from fewest to most generators
active_gen = sum(combinations(:,inc)>0,2);
[n_gen,sorted_i] = sort(active_gen);
combinations = combinations(sorted_i,:);

bal_field = fieldnames(qp.Organize.Balance);
n_f = length(bal_field);
for f = 1:1:n_f
    abbrev = nn_abrev{nonzeros((1:length(nn_list))'.*strcmp(bal_field{f},nn_list))};
    gen_dem = 0;
    if isfield(net_demand,abbrev)
        gen_dem = net_demand.(abbrev);
    end
    req = qp.Organize.Balance.(bal_field{f}); %rows of Aeq associated with electric demand
    [min_prod,max_prod] = check_feas(gen,combinations,options,qp,req,first_profile,dt,abbrev,false);
    keep = ((max_prod>=gen_dem) & (min_prod<=gen_dem));%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
    combinations = combinations(keep,:);
    n_gen = n_gen(keep);
end
end%ends function build_cases

function [min_prod,max_prod] = check_feas(gen,combinations,options,qp,req,ec,dt,abbrev,ar)
%returns the minimum and maximum production on the network speciied by
%abbrev, for each combination of generators specified by combinations
nn_list = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nn_abrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
n = length(combinations(:,1));
min_prod = zeros(n,length(gen));
max_prod = zeros(n,length(gen));

for i = 1:1:length(gen)
    states = qp.Organize.States{i};
    if ~isempty(states) && any(any(qp.Aeq(req,states)~=0))
        row_i = req(1);
        if length(req)>1
            net = nonzeros((1:length(nn_list))'.*strcmp(abbrev,nn_abrev));
            row_i = req(gen(i).QPform.(nn_list{net}).subnetNode);%need only the row corresponding to the node that this generator is on
        end
        switch gen(i).Type
            case 'Utility'%if it is a utility
                if length(states) ==2
                    min_prod(:,i) = sum(qp.Aeq(row_i,states(2))*qp.ub(states(2)));%if there is sellback, take 2nd state
                    max_prod(:,i) = sum(qp.Aeq(row_i,states(1))*qp.ub(states(1)));
                else
                    min_prod(:,i) = sum(qp.Aeq(row_i,states).*qp.lb(states)');
                    max_prod(:,i) = sum(qp.Aeq(row_i,states).*qp.ub(states)');
                end
            case {'Hydro Storage';'Electric Storage';'Thermal Storage';'Hydrogen Storage';}%if it is storage
                s = states(1);
                min_prod(:,i) = -gen(i).QPform.Ramp.b(1);%sum(qp.Aeq(row_i,s)*qp.lb(s));
                max_prod(:,i) = gen(i).QPform.Ramp.b(2);%sum(qp.Aeq(row_i,s)*qp.ub(s));
                if ~isempty(ec)%bounds reduced by state of charge
                    max_prod(:,i) = min(ec(i).*qp.Aeq(row_i,s)/dt,max_prod(:,i));
                    chargingSpace = sum((gen(i).QPform.Stor.UsableSize-ec(i)).*qp.Aeq(row_i,s));
                    min_prod(:,i) = max(-chargingSpace/dt, min_prod(:,i));
                end
            case 'AC_DC'
                if ~ar
                    switch abbrev
                        case 'E'
                            row_DC = gen(i).QPform.DirectCurrent.subnetNode;%need only the row corresponding to the node that this generator is on
                            [min_prod_DC,max_prod_DC] = check_feas(gen,combinations,options,qp,row_DC,ec,dt,'DC',true);
                            min_prod(:,i) = -gen(i).QPform.A.ub;%capacity to convert AC to DC power
                            max_prod(:,i) = qp.Aeq(row_i,states(2))*max_prod_DC;%pull out the DC-AC conversion efficiency
                        case 'DC'
                            row_AC = gen(i).QPform.Electrical.subnetNode;%need only the row corresponding to the node that this generator is on
                            [min_prod_AC,max_prod_AC] = check_feas(gen,combinations,options,qp,row_AC,ec,dt,'E',true);
                            min_prod(:,i) = -gen(i).QPform.B.ub;%capacity to convert DC to AC Power
                            max_prod(:,i) = qp.Aeq(row_i,states(1))*max_prod_AC;%pull out the AC-DC conversion efficiency
                    end
                end
            case {'Electric Generator';'CHP Generator';'Heater';'Electrolyzer';'Hydrogen Generator';'Chiller';'Cooling Tower';}
                if qp.Aeq(row_i,states(1))>0
                    min_prod(:,i) = (combinations(:,i)>0).*sum(qp.Aeq(row_i,states).*qp.lb(states)');
                    max_prod(:,i) = (combinations(:,i)>0).*sum(qp.Aeq(row_i,states).*qp.ub(states)');
                else%if it consumes this demand (i.e. a chiller when examining electrical demand), then lower_bound is -1*ub
                    min_prod(:,i) = (combinations(:,i)>0).*sum(qp.Aeq(row_i,states).*qp.ub(states)');
                    max_prod(:,i) = (combinations(:,i)>0).*sum(qp.Aeq(row_i,states).*qp.lb(states)');
                end
        end
        if isfield(gen(i).QPform,'constDemand') && isfield(gen(i).QPform.constDemand,abbrev)
            max_prod(:,i) = max_prod(:,i) - (combinations(:,i)>0).*gen(i).QPform.constDemand.(abbrev);
            min_prod(:,i) = min_prod(:,i) - (combinations(:,i)>0).*gen(i).QPform.constDemand.(abbrev);
        end
    end
end
if (strcmp(abbrev,'H') && options.excessHeat) || (strcmp(abbrev,'C') && options.excessCool)
    min_prod = min_prod - inf;
end

min_prod = sum(min_prod,2);
max_prod= sum(max_prod,2);
end%ends function check_feas

function cost = cost_estimate(c_red,gen,cost,verified,marginal,stor_pen,stor_abbrev,combinations,k)
%estimate the best achievable cost of an additional generator by adding
%its cost per kW at peak efficiency and subtracting the current marginal
%cost * the new generators capacity at peak efficiency
on_test = nonzeros(combinations(k,:));
n_on = length(on_test);
check = nonzeros((1:length(cost))'.*(~verified & sum(combinations(:,on_test)>0,2)==n_on & (isnan(cost) | cost<cost(k))));
d_cost = adjust_pen(stor_pen,gen,stor_abbrev,c_red,marginal);
for c = 1:1:length(check)
    i = check(c);
    new_gen = nonzeros(combinations(i,:)-combinations(k,:));
    cost(i) = max(cost(i),cost(k)+sum(d_cost(new_gen)));%The worst best case scenario is closest to the actual best case scenario, hence using max
end
end%Ends function cost_estimate

function d_cost = adjust_pen(stor_pen,gen,stor_abbrev,c_red,marginal)
n_g = length(gen);
d_cost = zeros(n_g,1);
d_pen = zeros(n_g,1);
for j = 1:1:length(stor_pen(1,:))
    x = stor_pen(1,j);
    H = stor_pen(2,j);
    if x>0 
        for i = 1:1:n_g
        	if c_red.p(i)>0 && isfield(gen(i).QPform.output,stor_abbrev{j}) && gen(i).QPform.output.(stor_abbrev{j})(1,1)>0%penalty can only be reduced by adding generators (of the correct type) if it was dawing power from the storage
                x_new = max(0,x - gen(i).Size*gen(i).QPform.output.(stor_abbrev{j})(1,1));%best case reduction in storage use
                d_pen(i) = d_pen(i) - 0.5*x^2*H + 0.5*x_new^2*H;
            end
        end
    end
end 
for i = 1:1:n_g
    if c_red.p(i)>0
        d_cost(i) = d_pen(i) - marginal.(c_red.ab{i})*c_red.p(i) + c_red.p(i)*c_red.c(i)/c_red.eff(i);
    end
end
end%ends function adjust_pen

function c_red = best_eff(gen,scale_cost,marginal)
n_g = length(gen);
c_red.c = zeros(n_g,1);
c_red.ab = cell(n_g,1);
c_red.eff = zeros(n_g,1);
c_red.p = zeros(n_g,1);
for i = 1:1:n_g
    switch gen(i).Type
        case {'CHP Generator';'Electric Generator';'Hydrogen Generator'}
            if isfield(gen(i).Output,'Electricity')
                [c_red.eff(i),p_i] = max(gen(i).Output.Electricity);
                c_red.ab{i} = 'E';
            else
                [c_red.eff(i),p_i] = max(gen(i).Output.DirectCurrent);
                c_red.ab{i} = 'DC';
            end
            c_red.c(i) = scale_cost(i);
            c_red.p(i) = gen(i).Output.Capacity(p_i)*gen(i).Size;
        case 'Heater'
            [c_red.eff(i),p_i] = max(gen(i).Output.Heat);
            c_red.ab{i} = 'H';
            c_red.c(i) = scale_cost(i);
            c_red.p(i) = gen(i).Output.Capacity(p_i)*gen(i).Size;
        case 'Chiller'
            [c_red.eff(i),p_i] = max(gen(i).Output.Cooling);
            c_red.ab{i} = 'C';
            if strcmp(gen(i).Source,'Heat')
                c_red.c(i) = marginal.H;
            else
                c_red.c(i) = marginal.E;
            end
            c_red.p(i) = gen(i).Output.Capacity(p_i)*gen(i).Size;
        case 'Cooling Tower'

        case 'Electrolyzer'


    end
end
end%Ends function best_eff