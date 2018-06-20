function solution = sort_solution(x,qp)
[m,n] = size(qp.organize);
n_s = m-1;
n_g = length(qp.constCost(1,:));
n_b = length(qp.Organize.Building.r(1,:));
n_l = length(qp.Organize.Transmission);
n_h = nnz(qp.Organize.Hydro);
n_fl = length(qp.Organize.fluid_loop);
solution.Dispatch = zeros(n_s+1,n_g);
solution.hydroSOC = zeros(n_s,n_h);
solution.excessHeat = [];
solution.excessCool = [];
solution.LineFlows = zeros(n_s,n_l);
solution.LineLoss = zeros(n_s,n_l);
solution.Buildings.Heating = zeros(n_s,n_b);
solution.Buildings.Cooling = zeros(n_s,n_b);
solution.Buildings.Temperature = zeros(n_s,n_b);
solution.fluid_loop = zeros(n_s,n_fl);
for i = 1:1:n_g
    if isfield(qp,'Renewable') && any(qp.Renewable(:,i)~=0)
        solution.Dispatch(2:end,i) = qp.Renewable(:,i);
    else
        out_vs_state = qp.Organize.Out_vs_State{1,i};%linear maping between state and output
        for t = 1:1:n_s+1
            if ~isempty(qp.organize{t,i})
                states = qp.organize{t,i};
                p = 0;
                for j = 1:1:length(states)
                    p = p + out_vs_state(j)*x(states(j)); %record this combination of outputs (SOC for storage)
                end
                if abs(p)>2e-4%added this to avoid rounding errors in optimization. When generator is locked off, the UB is set to 1e-5
                    solution.Dispatch(t,i) = p;
                end
            end
        end
    end
end
for i = 1:1:length(qp.Organize.Hydro)%Get SOC of each generator into a matrix for all time steps
    for t = 1:1:n_s
        solution.hydroSOC(t,i) = x(qp.organize{t+1,qp.Organize.Hydro(i)}+1,1)+qp.Organize.hydroSOCoffset(i);
    end
end
for i = 1:1:n_l
    for t = 1:1:n_s
        solution.LineFlows(t,i) = sum(x(qp.organize{t+1,i+n_g}));%power transfer or water flow rate
        if length(qp.Organize.States{i+n_g})>1
            solution.LineLoss(t,i) = sum(x(qp.organize{t+1,i+n_g}+1)); %down (positive) lines
            solution.LineLoss(t,i) = solution.LineLoss(t,i) + sum(x(qp.organize{t+1,i+n_g}+2)); %up (negative) lines
        end
    end
end
for i = 1:1:n_b
    for t = 1:1:n_s
        solution.Buildings.Temperature(t,i) = x(qp.organize{t,n_g+n_l+i},1);
        solution.Buildings.Heating(t,i) = x(qp.organize{t,n_g+n_l+i}+1,1) - qp.Organize.Building.H_Offset(t,i);
        solution.Buildings.Cooling(t,i) = x(qp.organize{t,n_g+n_l+i}+2,1) - qp.Organize.Building.C_Offset(t,i);
    end
end
for i = 1:1:n_fl
    for t = 1:1:n_s
        solution.fluid_loop(t,i) = x(qp.organize{t,i+n_g+n_l+n_b},1);
    end
end

solution.Dispatch(abs(solution.Dispatch)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors
%pull out any dumped heat
if ~isempty(qp.Organize.HeatVented)
    solution.excessHeat = zeros(n_s,length(qp.Organize.HeatVented(1,:)));
    for i = 1:1:length(qp.Organize.HeatVented(1,:))
        for t = 1:1:n_s
            if qp.Organize.HeatVented(t,i)>0
                solution.excessHeat(t,i) = x(qp.Organize.HeatVented(t,i));
            end
        end
    end
end
%pull out any dumped cooling
if ~isempty(qp.Organize.CoolVented)
    solution.excessCool = zeros(n_s,length(qp.Organize.CoolVented(1,:)));
    for i = 1:1:length(qp.Organize.CoolVented(1,:))
        for t = 1:1:n_s
            if qp.Organize.CoolVented(t,i)>0
                solution.excessCool(t,i) = x(qp.Organize.CoolVented(t,i));
            end
        end
    end
end
end%Ends function sort_solution