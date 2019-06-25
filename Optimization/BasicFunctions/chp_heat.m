function heat_out = chp_heat(gen,output)
n_s = length(output);
heat_out  = - gen.constDemand.H*(output>0);  
states = gen.states(1:nnz(~cellfun('isempty',gen.states(:,end))),end);
for t = 1:1:n_s
    net_out = 0;
    j = 1;
    while j<=length(states) && net_out<output(t)
        seg = min(output(t) - net_out,gen.(states{j}).ub(2));
        net_out = net_out + seg;
        heat_out(t) = heat_out(t) + seg*gen.output.H(j,2);
        j = j+1;
    end
end
end%Ends function chp_heat