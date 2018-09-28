function qp = disable_generators(qp,locked)
% Organize is the record of which indices are associated with each generator
% Locked is the specific combination of generators being tested
[m,n] = size(qp.organize);
n_s = max(1,m-1);
n_g = length(qp.Organize.Dispatchable);
% n_fl = length(organize.fluid_loop);
if isempty(locked)
    locked = true(m,n_g);
end
for i = 1:1:n_g
    if ~qp.Organize.Enabled(i)
        locked(:,i)=false;
    end
end

%add in constant electric loads from chillers/pumps in FitB
if isfield(qp,'constDemand') 
    outs = fieldnames(qp.constDemand);
    for s = 1:1:length(outs)
        if any(qp.constDemand.(outs{s}).req>0)
            for node= 1:1:qp.network.(outs{s}).nodes
                index = (node-1)*n_s+1:node*n_s;
                qp.constDemand.(outs{s}).load(index,:) = qp.constDemand.(outs{s}).load(index,:).*locked(2:end,:);
            end
            qp.beq(qp.constDemand.(outs{s}).req) = qp.beq(qp.constDemand.(outs{s}).req) + sum(qp.constDemand.(outs{s}).load,2);
        end
    end
end

%as you remove states and constraints the same variable at different time
%step will no longer be t1states away. Need to keep track of contraints for
%each time step seprately
f_name = fieldnames(qp.Organize.Balance);
for f = 1:1:length(f_name)
    n_bal = length(qp.Organize.Balance.(f_name{f}));
    qp.Organize.Balance.(f_name{f}) = ones(n_s,1)*qp.Organize.Balance.(f_name{f})' + (0:1:n_s-1)'*ones(1,n_bal)*qp.Organize.t1Balances;
end
if ~isempty(qp.Organize.SpinReserveStates)
    sr = true;
    qp.Organize.SpinReserve = [qp.Organize.SpinReserve;zeros(n_s-1,n_g+2);];
else
    sr = false;
end
qp.Organize.Ramp_up = zeros(n_s,n_g);
qp.Organize.Ramp_up(1,:) = qp.Organize.Ramping;
qp.Organize.Ramp_down = qp.Organize.Ramp_up+1;
a = qp.Organize.Equalities;
qp.Organize.Equalities = zeros(n_s,n_g,2);
qp.Organize.Equalities(1,:,:) = a;
b = qp.Organize.Inequalities;
qp.Organize.Inequalities = zeros(n_s,n_g,2);
qp.Organize.Inequalities(1,:,:) = b;
% c = qp.Organize.fluid_loop;
% qp.Organize.fluid_loop = zeros(n_s,n_fl);
% qp.Organize.fluid_loop(1,:) = c';
qp.Organize.SpinRow = [qp.Organize.SpinRow;zeros(n_s-1,n_g)];
n_b = length(qp.Organize.Building.r);
qp.Organize.Building.r = [qp.Organize.Building.r;zeros(n_s-1,n_b);];
if ~isempty(qp.Organize.HeatVented)
    qp.Organize.HeatVented = [qp.Organize.HeatVented;zeros(n_s-1,length(qp.Organize.HeatVented))];
end
if ~isempty(qp.Organize.CoolVented)
    qp.Organize.CoolVented = [qp.Organize.CoolVented;zeros(n_s-1,length(qp.Organize.CoolVented))];
end
for t = 2:1:n_s
    qp.Organize.Ramp_up(t,:) = qp.Organize.Ramp_up(t-1,:) + (qp.Organize.Ramp_up(1,:)>0)*qp.Organize.t1ineq;
    qp.Organize.Equalities(t,:,1) = qp.Organize.Equalities(t-1,:,1) + (qp.Organize.Equalities(1,:,1)>0)*qp.Organize.t1Balances;
    qp.Organize.Equalities(t,:,2) = qp.Organize.Equalities(t-1,:,2);
%     qp.Organize.fluid_loop(t,:) = qp.Organize.fluid_loop(t-1,:) + (qp.Organize.fluid_loop(1,:)>0)*qp.Organize.t1Balances;
    qp.Organize.Inequalities(t,:,1) = qp.Organize.Inequalities(t-1,:,1) + (qp.Organize.Inequalities(1,:,1)>0)*qp.Organize.t1ineq;
    qp.Organize.Inequalities(t,:,2) = qp.Organize.Inequalities(t-1,:,2);
    qp.Organize.SpinRow(t,:) = qp.Organize.SpinRow(t-1,:) + (qp.Organize.SpinRow(1,:)>0)*qp.Organize.t1ineq;
    qp.Organize.Building.r(t,:) = qp.Organize.Building.r(t-1,:) + (qp.Organize.Building.r(1,:)>0)*qp.Organize.t1ineq;
    if sr
        qp.Organize.SpinReserve(t,:) = qp.Organize.SpinReserve(t-1,:) + (qp.Organize.SpinReserve(1,:)>0)*qp.Organize.t1ineq;
        qp.Organize.SpinReserveStates(t,:) = qp.Organize.SpinReserveStates(t-1,:) + (qp.Organize.SpinReserveStates(1,:)>0)*qp.Organize.t1States;
    end
    if ~isempty(qp.Organize.HeatVented)
        qp.Organize.HeatVented(t,:) = qp.Organize.HeatVented(t-1,:) + (qp.Organize.HeatVented(1,:)>0)*qp.Organize.t1States;
    end
    if ~isempty(qp.Organize.CoolVented)
        qp.Organize.CoolVented(t,:) = qp.Organize.CoolVented(t-1,:) + (qp.Organize.CoolVented(1,:)>0)*qp.Organize.t1States;
    end
end
qp.Organize.Ramp_down = qp.Organize.Ramp_up+1*(qp.Organize.Ramp_up>0);

if any(any(~locked))
    [r,~] = size(qp.A);
    [req,x_l] = size(qp.Aeq);
    r_keep = true(r,1);
    req_keep = true(req,1);
    x_keep = true(x_l,1);
    rm_states = zeros(m,n);
    %%deal with initial condition
    for i = 1:1:n
        if i<=n_g &&  ~locked(1,i) 
            if ~locked(2,i)
                x_keep(qp.organize{1,i}) = false; %removes states associated with this generator.
                qp.organize{1,i} = [];
                rm_states(1,i:end) = rm_states(1,i:end) + 1;
                rm_states(2:end,:) = rm_states(2:end,:) + 1;
                req_keep(qp.Organize.IC(i)) = false;
                qp.Organize.IC(i) = 0;
            else %if it turns on right away keep IC = 0 constraint, but change ramping constraints
                if qp.Organize.Ramp_up(1,i)>0
                    s = qp.organize{1,i};
                    s2 = qp.organize{2,i}(1);
                    r_keep(qp.Organize.Ramp_down(1,i)) = false;%eliminate ramp down constraint
                    qp.b(qp.Organize.Ramp_up(1,i)) = max(qp.b(qp.Organize.Ramp_up(1,i)),-1.001*qp.A(qp.Organize.Ramp_up(1,i),s)*qp.lb(s2));% ramp up constraint
                    qp.lb(s) = 0;
                    qp.ub(s) = 0;
                end
                qp.beq(qp.Organize.IC(i)) = 0;
                qp.organize{1,i} = qp.organize{1,i} - rm_states(1,i);
                qp.Organize.IC(i) = qp.Organize.IC(i) - rm_states(1,i);
            end
        elseif ~isempty(qp.organize{1,i})
            qp.organize{1,i} = qp.organize{1,i} - rm_states(1,i);
            qp.Organize.IC(i) = qp.Organize.IC(i) - rm_states(1,i);
        end
    end  
    for t = 1:1:n_s
        for i = 1:1:n
            if i<=n_g && ~locked(t+1,i)
                states = qp.Organize.States{i};
                s_rm = [];
                if (t == n_s && ~locked(t,i)) || (~locked(t,i) && ~locked(t+2,i)) || qp.Organize.Ramping(1,i)==0%completely remove all generator states at this time step
                    s_rm = qp.organize{t+1,i}(1) + (0:(length(states)-1));
                    qp.organize{t+1,i} = [];
                else %keep the first state for ramping constraint to zero
                    if length(qp.organize{t+1,i})>1
                        s_rm = qp.organize{t+1,i}(2) + (0:(length(states)-2));
                    end
                    qp.organize{t+1,i} = qp.organize{t+1,i}(1) - rm_states(t+1,i); %;
                end
                if ~isempty(qp.Organize.SpinReserveStates) && qp.Organize.SpinReserveStates(t,i)>0%eliminate the spinning reserve state for the individual generator
                    s_rm = [s_rm qp.Organize.SpinReserveStates(t,i)];
                    qp.Organize.SpinReserveStates(t,i) = 0;
                end
                x_keep(s_rm) = false; %removes states associated with this generator.
                %keep track of the cumulative states removed
                rm_states(t+1,i:end) = rm_states(t+1,i:end) + length(s_rm);
                if t<n_s
                    rm_states(t+2:end,:) = rm_states(t+2:end,:) + length(s_rm);
                end
                if ~isempty(qp.Organize.Equalities) && qp.Organize.Equalities(t,i,1)>0%remove any associated equality constraints
                    eq = qp.Organize.Equalities(t,i,:);
                    req_keep(eq(1):eq(1)+eq(2)-1) = false;%equality constraints
                    qp.Organize.Equalities(t,i,1:2) = [0,0];
                end
                if qp.Organize.Ramp_up(t,i)>0 %remove the ramping constraint
                    if isempty(qp.organize{t+1,i})
                        r_keep(qp.Organize.Ramp_up(t,i)) = false;% ramp up constraint
                        r_keep(qp.Organize.Ramp_down(t,i)) = false;% ramp down constraint
                        if t<n_s
                            r_keep(qp.Organize.Ramp_up(t+1,i)) = false;% ramp up constraint
                            r_keep(qp.Organize.Ramp_down(t+1,i)) = false;% ramp down constraint
                        end
                    else %simplify the ramping constraint when ramping to/from zero
                        if i == 1
                            s = qp.organize{t+1,i}+rm_states(t,end);%state that is simply set equal to zero
                        else
                            s = qp.organize{t+1,i}+rm_states(t+1,i-1);%state that is simply set equal to zero
                        end
                        r_up = qp.Organize.Ramp_up(t,i);
                        r_down = qp.Organize.Ramp_down(t,i);
                        r_keep(r_up) = false;% ramp up constraint from previous condition
                        if ~locked(t,i)
                            r_keep(r_down) = false;% ramp down constraint
                        else
                            qp.b(r_down) = max(qp.b(r_down),-1.001*qp.A(r_down,s)*qp.lb(s));% ramp down constraint
                        end
                        if t == n_s
                            %no ramping constraint to next step to remove
                        elseif ~locked(t+2,i) %can remove subsequent ramp up & down constraint
                            r_keep(qp.Organize.Ramp_up(t+1,i)) = false;% ramp up constraint to next time
                            r_keep(qp.Organize.Ramp_down(t+1,i)) = false;% ramp down constraint to next time
                        else
                            r_keep(qp.Organize.Ramp_down(t+1,i)) = false;% ramp down constraint to next time
                            r_up2 = qp.Organize.Ramp_up(t+1,i);
                            qp.b(r_up2) = max(qp.b(r_up2),-1.001*qp.A(r_up2,s)*qp.lb(s));% ramp up constraint
                        end
                        qp.lb(s) = 0;
                        qp.ub(s) = 0;
                    end
                end
                if qp.Organize.SpinRow(t,i)>0
                    r_keep(qp.Organize.SpinRow(t,i)) = false;% spinning reserve inequality constraint 1 & 2
                    r_keep(qp.Organize.SpinRow(t,i)+1) = false;% spinning reserve inequality constraint 1 & 2
                end
                if ~isempty(qp.Organize.Inequalities) && qp.Organize.Inequalities(t,i,1)>0
                    ineq = qp.Organize.Inequalities(t,i,:);
                    r_keep(ineq(1):ineq(1)+ineq(2)-1) = false;%inequality constraints
                end
            else
                if ~isempty(qp.organize{t+1,i})%don't remove, but change state #'s
                    qp.organize{t+1,i} = qp.organize{t+1,i} - rm_states(t+1,i); %
                end
            end
        end
        %%update all other tracking variables that are used in sort solution
        if ~isempty(qp.Organize.HeatVented)
            for i = 1:1:length(qp.Organize.HeatVented(1,:))
                qp.Organize.HeatVented(t,i) = qp.Organize.HeatVented(t,i) - rm_states(t+1,n_g); 
            end
        end
        if ~isempty(qp.Organize.CoolVented)
            for i = 1:1:length(qp.Organize.CoolVented(1,:))
                qp.Organize.CoolVented(t,i) = qp.Organize.CoolVented(t,i) - rm_states(t+1,n_g); 
            end
        end
        if sr
            for i = 1:1:length(qp.Organize.SpinReserveStates(1,:))
                if qp.Organize.SpinReserveStates(t,i)>0
                    qp.Organize.SpinReserveStates(t,i) = qp.Organize.SpinReserveStates(t,i) - rm_states(t+1,min(i,n_g));
                end
            end
        end
        for f = 1:1:length(f_name)
            qp.Organize.Balance.(f_name{f})(t,:) = qp.Organize.Balance.(f_name{f})(t,:) - nnz(~req_keep(1:qp.Organize.Balance.(f_name{f})(t,1)));
        end
        %the following tracks the constraint rows, but is probably not necessary
        for i = 1:1:n_g
            if qp.Organize.Ramp_up(t,i)>0 
                if r_keep(qp.Organize.Ramp_up(t,i))
                    qp.Organize.Ramp_up(t,i) = qp.Organize.Ramp_up(t,i) - nnz(~r_keep(1:qp.Organize.Ramp_up(t,i)));
                else
                    qp.Organize.Ramp_up(t,i) = 0;
                end
            end
            if qp.Organize.Ramp_down(t,i)>0 
                if r_keep(qp.Organize.Ramp_down(t,i))
                    qp.Organize.Ramp_down(t,i) = qp.Organize.Ramp_down(t,i) - nnz(~r_keep(1:qp.Organize.Ramp_down(t,i)));
                else
                    qp.Organize.Ramp_down(t,i) = 0;
                end
            end
            if ~isempty(qp.Organize.Equalities) && qp.Organize.Equalities(t,i,1)>0
                qp.Organize.Equalities(t,i,1) = qp.Organize.Equalities(t,i,1) - nnz(~req_keep(1:qp.Organize.Equalities(t,i,1)));
            end
            if qp.Organize.SpinRow(t,i)>0
                qp.Organize.SpinRow(t,i) = qp.Organize.SpinRow(t,i) - nnz(~r_keep(1:qp.Organize.SpinRow(t,i)));
            end
            if ~isempty(qp.Organize.Inequalities) && qp.Organize.Inequalities(t,i)>0
                qp.Organize.Inequalities(t,i) = qp.Organize.Inequalities(t,i) - nnz(~r_keep(1:qp.Organize.Inequalities(t,i)));
            end
        end
%         for i = 1:1:n_fl
%             if ~isempty(qp.Organize.fluid_loop) && qp.Organize.fluid_loop(t,i)>0
%                 qp.Organize.fluid_loop(t,i) = qp.Organize.fluid_loop(t,i) - nnz(~req_keep(1:qp.Organize.fluid_loop(t,i)));
%             end
%         end
    end
    %eliminate states and constraints so that problem is in reduced form
    qp.H = qp.H(x_keep,x_keep);
    qp.f = qp.f(x_keep);
    qp.Aeq = qp.Aeq(req_keep,x_keep);
    qp.beq = qp.beq(req_keep);
    if r>0 && any(r_keep)
        qp.A = qp.A(r_keep,x_keep);
        qp.b = qp.b(r_keep);
    else
        qp.A = [];
        qp.b = [];
    end
    qp.lb = qp.lb(x_keep);
    qp.ub = qp.ub(x_keep);
end
end%Ends function disable_generators