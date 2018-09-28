function solution = single_optimization_NREL(gen,buildings,options,test_data,date1)
n_g = length(gen);%skip initialization
for i=1:1:n_g
    if ismember(gen(i).Type,{'Electric Storage';'Thermal Storage';})
        gen(i).CurrentState(1) = 0.5*gen(i).Size*(gen(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
    else 
        gen(i).CurrentState(1) = 0;
    end
end

test_data = update_test_data(test_data,[],[],options);
date = date1+[0;build_time_vector(options)/24];%linspace(Date,DateEnd)';would need to re-do optimization matrices for this time vector
[forecast,gen,~] = update_forecast(gen,buildings,[],[],options,date(2:end),test_data,[]);
scaleCost = update_cost(date(2:end),gen); %% All feedstock costs were assumed to be 1 when building matrices 
solution = solver_nrel(gen,options,22,forecast,scaleCost);
solution.Timestamp = date;
end %Ends function single_optimization_NREL