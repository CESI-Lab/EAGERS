function [max_out,constraint,spare_gen_cumulative] = gen_limit(gen,gen_output,binary,dt)
%% Find Max output of each generator at each time
n_g = length(gen);
n_s = length(dt);
ramp_up = zeros(n_s,n_g);
ramp_down = zeros(n_s,n_g);
upper_bound = zeros(1,n_g);
lower_bound = zeros(1,n_g);
constraint = zeros(n_s,n_g);
for i = 1:1:n_g
    if isfield(gen(i).QPform,'Ramp')
        ramp_up(:,i) = gen(i).QPform.Ramp.b(1)*dt;
        ramp_down(:,i) = gen(i).QPform.Ramp.b(2)*dt;
        starts = nonzeros((1:n_s)'.*((binary(2:end,i)-binary(1:n_s,i))>0));
        stops = nonzeros((1:n_s)'.*(binary(1:n_s,i)-(binary(2:end,i))>0));
        lower_bound(i) = gen(i).QPform.(gen(i).QPform.states{1,end}).lb(end);
        ramp_up(starts,i) = max(ramp_up(starts,i),lower_bound(i)); %if the ramp rate is less than the lb, increase ramp rate at moment of startup
        ramp_down(stops,i) = max(ramp_down(stops,i),lower_bound(i)); %if the ramp rate is less than the lb, increase ramp rate at moment of shutdown
        upper_bound(i) = gen(i).Size;
    end
end
s= {'E';'H';'C';'CW';'Hy'};
s2 = {{'CHP Generator';'Electric Generator';};{'Heater'};{'Chiller';};{'Chiller';'CoolingTower'};{'Electrolyzer';'Hydrogen Generator';}};
for k = 1:1:length(s)
    spare_gen = zeros(n_s,1);
    spare_gen_cumulative.(s{k}) = zeros(n_s,1);
    max_out.(s{k}) = zeros(n_s+1,n_g);
    inc = false(n_g,1);
    type = cell(n_g,1);
    for i = 1:1:n_g
        type{i} = gen(i).Type;
        if any(strcmp(gen(i).Type,s2{k}))
            inc(i) = true;
        end
    end
    for i = 1:1:n_g
        if inc(i)
            max_out.(s{k})(1,i) = gen_output(1,i);
            start = inf; %constrained by initial condition (can't do anything)
            for t = 1:1:n_s
                if binary(t+1,i)
                    if~binary(t,i)
                        start = t;%just turned on
                    end
                    if upper_bound(i)<=(max_out.(s{k})(t,i) + ramp_up(t,i))
                        max_out.(s{k})(t+1,i) = upper_bound(i);
                    else
                        max_out.(s{k})(t+1,i) = (max_out.(s{k})(t,i) + ramp_up(t,i));
                        if ~isinf(start) && start>1
                            constraint(t,i) = start-1;
                        else
                            constraint(t,i) = inf;
                        end
                    end
                end
            end
            for t = n_s:-1:2
                if max_out.(s{k})(t,i) - ramp_down(t,i) > max_out.(s{k})(t+1,i) %constrained by shutdown
                    max_out.(s{k})(t,i) = min(upper_bound(i),max_out.(s{k})(t+1,i) + ramp_down(t,i));
                    if constraint(t-1,i) == 0
                        constraint(t-1,i) = n_s +1 - max(((n_s-t+1):-1:1)'.*(max_out.(s{k})(t+1:n_s+1,i)==0));
                    end
                end
            end
            if strcmp(gen(i).Type,'Chiller') && strcmp(gen(i).Source,'Heat')%assume net heat production remains the same, so do not count spare absorption chiller capacity
                max_out.(s{k})(2:end,i) = gen_output(2:end,i);%don't count towards spare capacity
            end
            spare_gen = spare_gen + (max_out.(s{k})(2:end,i) - gen_output(2:end,i));
        end
        if strcmp(s{k},'H') && strcmp(gen(i).Type,'CHP Generator')
            heat_output = chp_heat(gen(i).QPform,gen_output(:,i));%Convert GenOutput to heatOutput
            max_out.H(:,i) = max_out.H(:,i) + chp_heat(gen(i).QPform,max_out.E(:,i));%convert MaxOut.E to MaxOut.H
            spare_gen = spare_gen + (max_out.(s{k})(2:end,i) - heat_output(2:end));
        end
        if strcmp(gen(i).Type,'Utility') && isfield(gen(i).QPform.output,s{k})
            max_out.(s{k})(2:end,i) = gen(i).QPform.X.ub; gen_output(2:end,i);
            spare_gen = spare_gen + gen(i).QPform.X.ub - gen_output(2:end,i);
        end
    end    
    for t = 1:1:n_s
        spare_gen_cumulative.(s{k})(t) = sum(spare_gen(1:t).*dt(1:t));
    end
    %%allow for spare DC generation to supply AC power needs & vice versa
    if strcmp(s{k},'E') && any(strcmp('AC_DC',type))
        acdc = nonzeros((1:1:n_g)'.*strcmp('AC_DC',type));
        max_out.DC = zeros(n_s+1,n_g);
        spare_gen_cumulative.DC = spare_gen_cumulative.E;
        for i = 1:1:n_g
            if inc(i) 
                if isfield(gen(i).QPform.output,'DC')
                    max_out.DC(:,i) = max_out.E(:,i);
                    max_out.E(:,acdc) = max_out.E(:,acdc) + max_out.E(:,i)*gen(acdc).VariableStruct.DC_to_AC_eff;
                else
                    max_out.DC(:,i) = max_out.E(:,i)*gen(acdc).VariableStruct.AC_to_DC_eff;
                end
            end
        end
    end
end
end%ends function gen_limit