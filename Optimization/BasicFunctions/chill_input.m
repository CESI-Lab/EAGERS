function input = chill_input(gen,output)
n_s = length(output);
input = zeros(n_s,1);
if strcmp(gen.Source,'Heat')
    coef = -gen.QPform.output.H;
    if isfield(gen.QPform.constDemand,'H')
        input = (output>1e-3)*gen.QPform.constDemand.H;
    end
else
    coef = -gen.QPform.output.E;
    if isfield(gen.QPform.constDemand,'E')
        input = (output>1e-3)*gen.QPform.constDemand.E;
    end
end  
n_state = length(gen.QPform.states(:,1));
for t = 1:1:n_s
    if output(t)>1e-3
        net_out = 0;
        i = 1;
        while i<= n_state && net_out<output(t)
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
end%ends function chill_input