function input = cool_tower_input(gen,output)
n_s = length(output);
input = zeros(n_s,1);
coef = -gen.QPform.output.E;
for t = 1:1:n_s
    if output(t)>1e-3
        i = 1;
        net_out = 0;
        if isfield(gen.QPform.constDemand,'E')
            input(t) = gen.QPform.constDemand.E;
        end
        while net_out<output(t)
            state = gen.QPform.states{i,end};
            seg = gen.QPform.(state).ub(end);
            seg = min(seg,output(t)-net_out);
            net_out = net_out + seg;
            if length(coef(1,:))>1
                input(t) = input(t) + seg*coef(i,2);
            else
                input(t) = input(t) + seg*coef;
            end
            i = i+1;
        end  
    end
end
end%ends function cool_tower_input