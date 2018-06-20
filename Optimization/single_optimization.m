function [solution,forecast,gen,buildings,fluid_loop] = single_optimization(date,last_solution,gen,buildings,fluid_loop,options,dispatch,subnet,op_mat_a,op_mat_b,one_step,hist_prof,wy_forecast,prev_data,now_data,future_data)
% This function runs a single optimization over the horizon of the dispatch optimization
% date is a vector of timestamps representing the points in time to be forecast and optimized
% last_solution is optional. If it is not empty, [], it is a structure of the solution 1 step prior
% It returns solution, a structure organizing the generator setpoints, line
% transfer flow, excess heat rejected, and hydro state of charge, and it
% returns the forecast that was optimized.
% If the Optimization has not been run previously it loads the optimization
% matrices and automatically determines an initial condition. If it has
% been run previously it uses the end point of a prior solution that is
% closest to the first date in ForecastTime as the initial condition. 
% This function then generates a forecast of the demands and optimizes
% according to that forecast.
if ~isempty(dispatch) && isfield(dispatch,'Timestamp') && any(dispatch.Timestamp==(date(1) - options.Resolution/24))
    t = max(linspace(1,length(dispatch.Timestamp))'.*(dispatch.Timestamp<date(1))); %index preceeding current step
    for i = 1:1:length(gen)
        gen(i).CurrentState = dispatch.GeneratorState(t,:);
    end
else
    [data_t0,gen,~] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(1),hist_prof,prev_data,now_data,future_data);
    [gen,fluid_loop] = automatic_ic(gen,buildings,fluid_loop,subnet,date(1),one_step,options,data_t0);% set the initial conditions
end
[forecast,gen,buildings] = update_forecast(gen,buildings,fluid_loop,subnet,options,date(2:end),hist_prof,prev_data,now_data,future_data);
forecast.wy_forecast = wy_forecast;
[solution,~] = dispatch_loop(gen,buildings,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,last_solution);
end%ends Function single_optimization
