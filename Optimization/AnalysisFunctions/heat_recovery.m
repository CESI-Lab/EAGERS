function heat_recovery = heat_recovery(gen,dispatch)
n_g = length(gen);
[n_s,~] = size(dispatch);
heat_recovery = zeros(n_s,1);
for i = 1:1:n_g
    if strcmp(gen(i).Type,'CHP Generator') || (strcmp(gen(i).Type,'Chiller') && strcmp(gen(i).Source,'Heat'))
        gen_i = gen(i).QPform;
        states = gen_i.states(1:nnz(~cellfun('isempty',gen_i.states(:,end))),end);
        for t = 1:1:n_s
            pow = 0;
            if dispatch(t,i)>0 && isfield(gen_i,'constDemand') && isfield(gen_i.constDemand,'H')
                heat_recovery(t) = heat_recovery(t) -gen_i.constDemand.H;
            end
            h_ratio = gen_i.output.H(:,end);
            for j = 1:1:length(states)
                if length(h_ratio)>1
                    h_r = h_ratio(j);
                else
                    h_r = h_ratio;
                end
                heat_recovery(t) = heat_recovery(t) + min(dispatch(t,i)-pow,gen_i.(states{j}).ub(end))*h_r;
                pow = pow+min(dispatch(t,i)-pow,gen_i.(states{j}).ub(end));
            end
        end
    end
end
heat_recovery = heat_recovery/293.07;%convert to mmBTU
end%ends function heat_recovery