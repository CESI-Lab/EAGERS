function [solution,lb_relax] = cqp_method(gen,qp_0,first_disp,date)
n_s = length(date)-1;
n_g = length(gen);
locked = true(n_s+1,n_g);
for i = 1:1:n_g
    if ~gen(i).Enabled
        locked(:,i) = 0;
    end
end
lower_bound = zeros(1,n_g);
upper_bound = zeros(1,n_g);
dx = zeros(n_s,n_g);
dt = (date(2:end) - date(1:end-1))*24;
for i = 1:1:n_g
    if qp_0.Organize.Dispatchable(i) ==1
        states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
        for j = 1:1:length(states)
            lower_bound(i) = lower_bound(i) + gen(i).QPform.(states{j}).lb(2);
            upper_bound(i) = upper_bound(i) + gen(i).QPform.(states{j}).ub(end);
        end
        locked(first_disp(:,i)<=lower_bound(i),i) = false;%Default to off when initial dispatch is below LB
    end
    if isfield(gen(i).VariableStruct,'dX_dt')
        dx(:,i) = dt*gen(i).VariableStruct.dX_dt;
    end
end
%% make sure it can shut down in time from initial condition
for i = 1:1:n_g
    if qp_0.Organize.Dispatchable(i) ==1
        if first_disp(1,i)>0 && ~all(locked(:,i))
            r = qp_0.Organize.Ramping(i)+1;
            d = first_disp(1,i);
            t = 1;
            while d>0
                d = d - qp_0.b(r);
                if d>0 && ~locked(t+1,i)
                    locked(t+1,i) = true;
                end
                t = t+1;
                r = r+qp_0.Organize.t1ineq;
            end
        end
    end
end

feasible = 0;
attempt = 0;
lb_relax = 1;

while feasible ~= 1 && attempt <8
    %attempt: integer value describing the number of attempts before reaching
    %feasibility, this determines how close components must be to their lower
    %bound from below to be considered online
    %n represents the percent of lower bounds 
    %on your first try, just use the locked matrix given, then do unit
    %commitment based on OptimalState>LB*percLB
    perc_lb = [0.9, 0.75, 0.5, 0.2, 0.1, 0, -1];

    if attempt>0%second try, lower limit for online threshold
        lb_relax = perc_lb(attempt);
        %only change label for unit commitment gens, and don't change the
        %label for initial conditions
        for i = 1:1:n_g
            if qp_0.Organize.Dispatchable(i) ==1
                locked(2:end,i) = (first_disp(2:end,i)>(lower_bound(i)*lb_relax));%Default to on unless offline in initial dispatch
            end
        end
    end
    qp = disable_generators(qp_0,locked);%Disable generators here
    [x,feasible] = call_solver(qp);
    attempt = attempt+1;
end


if feasible==1
    solution = sort_solution(x,qp);
    first_disp = solution.Dispatch;
    cost = sum(net_cost(gen,solution.Dispatch,date,'Dispatch'));
    for i = 1:1:n_g
        if qp_0.Organize.Dispatchable(i)
            [locked,cost,solution,starts] = rule_one(gen,locked,cost,solution,date,qp_0,dx,first_disp,i);
            [locked,cost,solution,stops] = rule_two(gen,locked,cost,solution,date,starts,qp_0,dx,first_disp,i);
            [locked,cost,solution] = rule_three(gen,locked,cost,solution,date,starts,stops,qp_0,dx,first_disp,i);
        end
    end
    [~,~,solution] = rule_four(gen,locked,cost,solution,date,qp_0,dx,upper_bound);
    %add something to see if FilterGenerators changed anything
else
    disp('cqp error: Cannot Find Feasible Dispatch');
end
solution.LBrelax = lb_relax;
end%ends function cqp_method

function [locked,cost,solution,starts] = rule_one(gen,locked,cost,solution,date,qp_0,dx,first_disp,i)
%% Rule 1: turn off for longer time at start if possible
n_s = length(dx(:,1));
index = (1:n_s)';
starts = nonzeros(index.*((locked(2:end,i)-locked(1:n_s,i))>0)); % if on at t = 3, start = 3, where IC is t=0
if ~locked(1,i) && ~isempty(starts)
    p = 0;
    ramp_up = dx(starts(1),i);
    while ramp_up<first_disp(starts(1)+1,i) && (starts(1)-p>0)
        ramp_up = ramp_up+dx(starts(1)-p,i);
        p = p+1;
    end
    if starts(1)-p>0
        l2 = locked;
        l2(1:(starts(1)-p+1),i) = false;
        qp = disable_generators(qp_0,l2);%Disable generators here
        [x,Feasible] = call_solver(qp);
        if Feasible == 1
            sol_new = sort_solution(x,qp);
            new_cost = sum(net_cost(gen,sol_new.Dispatch,date,'Dispatch'));
            if new_cost<cost
                locked = l2;
                cost = new_cost;
                solution = sol_new;
            end
        end
    end
    if length(starts)>1
        starts = starts(2:end);
    else
        starts =[];
    end
end
end%ends function rule_one

function [locked,cost,solution,stops] = rule_two(gen,locked,cost,solution,date,starts,qp_0,dx,first_disp,i)
%% Rule 2: If off for a long enough segment in first dispatch, try turning off for as much of that as possible given ramp rates
n_s = length(dx(:,1));
index = (1:n_s)';
stops = nonzeros(index.*((locked(1:n_s,i)-locked(2:end,i))>0)); % if off at t = 12, stop = 12, where IC is t=0
for k = 1:1:length(starts)
    if ~isempty(stops) && length(stops)>=k
        if sum(dx(stops(k):starts(k)-1,i))>(first_disp(stops(k),i) + first_disp(starts(k)+1,i)) %can ramp all the way down, and back up, and be off for 1 step
            l2 = locked;
            %find step when it can hit zero given setting at Disp(stops(k)
            n=1;
            ramp_down = dx(stops(k),i);
            while ramp_down<first_disp(stops(k),i)
                ramp_down = ramp_down+dx(stops(k)+n,i);
                n = n+1;
            end
            p = 1;
            ramp_up = dx(starts(k),i);
            while ramp_up<first_disp(starts(k)+1,i)
                ramp_up = ramp_up+dx(starts(k)-p,i);
                p = p+1;
            end
            l2((stops(k)+n):(starts(k)-p+1),i) = false;
            qp = disable_generators(qp_0,l2);%Disable generators here
            [x,feasible] = call_solver(qp);
            if feasible == 1
                sol_new = sort_solution(x,qp);
                new_cost = sum(net_cost(gen,sol_new.Dispatch,date,'Dispatch'));
                if new_cost<cost
                    locked = l2;
                    cost = new_cost;
                    solution = sol_new;
                end
            end
        end
    end
end
end%ends function rule_two

function [locked,cost,solution] = rule_three(gen,locked,cost,solution,date,starts,stops,qp_0,dx,first_disp,i)
%% Rule 3: try turning off @ end
n_s = length(dx(:,1));
if length(stops)>length(starts)
    n=stops(end);
    ramp_down = dx(n,i);
    while ramp_down<first_disp(stops(end),i) && n<(n_s) && n>0
        ramp_down = ramp_down+dx(n,i);
        n = n-1;
    end
    if n<(n_s)
        l2 = locked;
        l2((n+1):n_s+1,i) = false;
        qp = disable_generators(qp_0,l2);%Disable generators here
        [x,feasible] = call_solver(qp);
        if feasible == 1
            sol_new = sort_solution(x,qp);
            new_cost = sum(net_cost(gen,sol_new.Dispatch,date,'Dispatch'));
            if new_cost<cost
                locked = l2;
                cost = new_cost;
                solution = sol_new;
            end
        end
    end
end
end%ends rule_three

function [locked,cost,solution] = rule_four(gen,locked,cost,solution,date,qp_0,dx,upper_bound)
%% Rule 4: if on for a short time and sufficient capacity in remaining active generators, turn gen off completely
n_s = length(dx(:,1));
n_g = length(gen);
index = (1:n_s)';
for i = 1:1:n_g
    starts = nonzeros(index.*((locked(2:end,i)-locked(1:n_s,i))>0)); % if on at t = 3, start = 3, where IC is t=0
    stops = nonzeros(index.*((locked(1:n_s,i)-locked(2:end,i))>0)); % if off at t = 12, stop = 12, where IC is t=0
    if qp_0.Organize.Dispatchable(i)
        if locked(1,i)
            if length(stops)>1
                stops = stops(2:end);
            else
                stops =[];
            end
        end
        for k = 1:1:length(stops)
            if sum(dx(starts(k):stops(k),i))<upper_bound(i) && sum(locked(:,i))<floor(n_s/4)%can only ramp to 1/2 power and less than 1/4 of horizon
                l2 = locked;
                l2(starts(k):(stops(k)+1),i)= false;
                qp = disable_generators(qp_0,l2);%Disable generators here
                [x,feasible] = call_solver(qp);
                if feasible == 1
                    sol_new = sort_solution(x,qp);
                    new_cost = sum(net_cost(gen,sol_new.Dispatch,date,'Dispatch'));
                    if new_cost<cost
                        locked = l2;
                        cost = new_cost;
                        solution = sol_new;
                    end
                end
            end
        end
    end
end
end%ends function rule_four
