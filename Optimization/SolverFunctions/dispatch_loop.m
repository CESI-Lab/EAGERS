function [solution,flag] = dispatch_loop(gen,building,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,last)
%DISPATCH_LOOP 
%
% Flag values:
%   0 -- Standard operation.
%   1 -- Initial dispatch (Fit A) infeasible.
%   2 -- Second dispatch (Fit B) infeasible.
%   3 -- ANN unit commitment infeasible.

%% Initialize flag.
flag = 0;

%% calculate optimal dispatch over the forecast horizon
dt = (date(2:end) - date(1:end-1))*24;
n_g= length(gen);
scale_cost = update_cost(date(2:end),gen); %% All feedstock costs were assumed to be 1 when building matrices 
%% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
predict_dispatch = zeros(length(date),n_g);
if isempty(last) || isempty(last.Dispatch)
    for i = 1:1:n_g
        predict_dispatch(:,i) = gen(i).CurrentState(1);
    end
else
    index = [(3:length(last.Dispatch(:,1)))';length(last.Dispatch(:,1));];
    for i = 1:1:n_g
        predict_dispatch(:,i) = [gen(i).CurrentState(1); last.Dispatch(index,i)];
    end
end
margin_cost = update_mc(gen,predict_dispatch,scale_cost,dt,[]);
if strcmp(options.solver,'Gurobi')
    qp = update_matrices(gen,building,fluid_loop,subnet,options,op_mat_b,date,scale_cost,margin_cost,forecast,[]);
    tic
    x = gurobi_opt(qp);
    solution = sort_solution(x,qp);
    tsim(1,1) = toc;
elseif strcmp(options.solver,'ANN')
    % Step 1 & 2: unit commitment, handled by a trained ANN
    tic
    %train the first time through
    training = isempty(last.Dispatch);
    locked = fireANN(last.Dispatch,forecast,scale_cost,training);
    qp = update_matrices(gen,building,fluid_loop,subnet,options,op_mat_b,date,scale_cost,margin_cost,forecast,[]);
    %make sure initial condition of locked is correct
    for i = 1:1:n_g
        if qp.Organize.Dispatchable(i) ==1
            locked(1,i) = gen(i).Status;
        end
    end
    tsim(1,2) = toc;
    % Step 3: complete optimization with fit B
    tic
    qp = disable_generators(qp,locked); % Disable generators here
    [x,feasible] = call_solver(qp);
    if feasible == 1
        solution = sort_solution(x,qp);
    else
        %was infeasible using ANN; default to mcQP
        disp(strcat('ANN unit commitment infeasible, defaulting to mcQP at timestep', num2str(date(1))))
        solution = dispatch_loop(gen,building,fluid_loop,subnet,op_mat_a,op_mat_b,one_step,options.SpinReserve,true,'quadprog',date,forecast,last);
    end
    tsim(1,3) = toc;
else
    qp_1 = update_matrices(gen,building,fluid_loop,subnet,options,op_mat_a,date,scale_cost,margin_cost,forecast,[]);
    %% Step 1: Determine initial dispatch
    tic
    qp = disable_generators(qp_1,[]);%Disable generators here
    [x,feasible] = call_solver(qp);
    if feasible == 1
        solution1 = sort_solution(x,qp);
    else
        %% add infinite sources for each energy balance and determine short-commings
        flag = 1;
        solution1 = fix_infeasible(qp);
    end
    tsim(1,1) = toc;
    if feasible == 1 && any(qp.Organize.Dispatchable) %might be some on/off combinations
        if options.MixedInteger
            tic
            %% Step 2: Dispatch step by step
            [optimal_state,flag] = dispatch_step(gen,building,fluid_loop,subnet,options,one_step,date(1),forecast,scale_cost,dt,solution1.Dispatch);
            flag = flag * 2;  % Convert: 1 -> 2.
            clear mex
            tsim(1,2) = toc;
        else
            optimal_state = solution1.Dispatch;
        end
        %% Step 3: 2nd complete optimization
        tic
        margin_cost = update_mc(gen,optimal_state,scale_cost,dt,[]);
        qp_2 = update_matrices(gen,building,fluid_loop,subnet,options,op_mat_b,date,scale_cost,margin_cost,forecast,[]); % Update fit B matrices.
        if flag == 2 || ~options.MixedInteger
            solution = cqp_method(gen,qp_2,solution1.Dispatch,date);
        else
            locked = verify_ramping(gen,subnet,qp_2,optimal_state,dt);
            qp = disable_generators(qp_2,locked); % Disable generators here.
            [x,feasible] = call_solver(qp); % This is the dispatch with fit B.
            if feasible == 1
                solution = sort_solution(x,qp);
            else
                flag = 3;
%                 disp('dispatch error: Cannot Find Feasible Dispatch');
                solution = solution1;
            end
        end
        tsim(1,3) = toc;
    else
        solution = solution1;
    end
end
solution.timers = tsim;
end % Ends function dispatch_loop