function [dispatch,line_loss,excess_heat,excess_cool,hydro_soc,temperature,heating,cooling] = sort_solution_step(x,qp)
n_g = length(qp.constCost);
n_b = length(qp.Organize.Building.r);
n_l = length(qp.Organize.Transmission);
n_fl = length(qp.Organize.fluid_loop);
n_h = nnz(qp.Organize.Hydro);
hydro_soc = zeros(1,n_h);
excess_heat = zeros(1,nnz(qp.Organize.HeatVented));
excess_cool = zeros(1,nnz(qp.Organize.CoolVented));
line_flows = zeros(1,n_l);
line_loss = zeros(1,n_l);
temperature = zeros(1,n_b);
fluid_loop = zeros(1,n_fl);
heating = zeros(1,n_b);
cooling = zeros(1,n_b);
gen_disp = zeros(1,n_g);
for i = 1:1:n_g
    if isfield(qp,'Renewable') && any(qp.Renewable(:,i)~=0)
        gen_disp(2:end,i) = qp.Renewable(1,i);
    else
        out_vs_state = qp.Organize.Out_vs_State{1,i};%linear maping between state and output
        if ~isempty(qp.organize{1,i})
            states = qp.organize{1,i};
            for j = 1:1:length(states)
                gen_disp(1,i) = gen_disp(1,i) + out_vs_state(j)*x(states(j)); %record this combination of outputs (SOC for storage)
            end
        end
    end
    
end
for i = 1:1:length(qp.Organize.Hydro)
    hydro_soc(1,i) = x(qp.organize{1,qp.Organize.Hydro(i)}+1,1);
end
gen_disp(abs(gen_disp)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors

for i = 1:1:n_l
    line_flows(1,i) = sum(x(qp.organize{1,i+n_g}));%power transfer
    if length(qp.Organize.States{i+n_g})>1
        line_loss(1,i) = sum(x(qp.organize{1,i+n_g}+1)); %down (positive) lines
        line_loss(1,i) = line_loss(1,i) + sum(x(qp.organize{1,i+n_g}+2)); %up (negative) lines
    end
end

for i = 1:1:n_b
    temperature(1,i) = x(qp.organize{1,i+n_g+n_l},1);
    heating(1,i) = x(qp.organize{1,i+n_g+n_l}+1,1) + qp.Organize.Building.H_Offset(1,i);
    cooling(1,i) = x(qp.organize{1,i+n_g+n_l}+2,1) + qp.Organize.Building.C_Offset(1,i);
end
for i = 1:1:n_fl
    fluid_loop(1,i) = x(qp.organize{1,i+n_g+n_l+n_b},1);
end
dispatch = [gen_disp,line_flows,temperature,fluid_loop];
%pull out any dumped heat
for i = 1:1:length(qp.Organize.HeatVented)
    if qp.Organize.HeatVented(i)>0
        excess_heat(1,i) = x(qp.Organize.HeatVented(i));
    end
end
%pull out any dumped cooling
for i = 1:1:length(qp.Organize.CoolVented)
    if qp.Organize.CoolVented(i)>0
        excess_cool(1,i) = x(qp.Organize.CoolVented(i));
    end
end
end%Ends function sort_solution_step