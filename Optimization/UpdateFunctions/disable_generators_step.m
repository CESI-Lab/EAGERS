function qp = disable_generators_step(qp,locked)
% Organize is the record of which indices are associated with each generator
% Locked is the specific combination of generators being tested
n = length(qp.organize);
n_g = length(qp.Organize.Dispatchable);
if isempty(locked)
    locked = true(1,n_g);
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
                qp.constDemand.(outs{s}).load(node,:) = qp.constDemand.(outs{s}).load(node,:).*locked;
            end
            qp.beq(qp.constDemand.(outs{s}).req) = qp.beq(qp.constDemand.(outs{s}).req) + sum(qp.constDemand.(outs{s}).load,2);
        end
    end
end
if any(~locked)
    [r,~] = size(qp.A);
    [req,x_l] = size(qp.Aeq);
    r_keep = true(r,1);
    req_keep = true(req,1);
    x_keep = true(x_l,1);
    rm_states = zeros(1,n);
    for i = 1:1:n
        if i<=n_g && ~locked(1,i)
            states = qp.Organize.States{i};
            s_rm = qp.organize{1,i}(1) + (0:(length(states)-1));
            qp.organize{1,i} = [];
            if ~isempty(qp.Organize.SpinReserveStates) && qp.Organize.SpinReserveStates(1,i)>0%eliminate the spinning reserve state for the individual generator
                s_rm = [s_rm,qp.Organize.SpinReserveStates(1,i)];
                qp.Organize.SpinReserveStates(1,i) = 0;
            end
            x_keep(s_rm) = false; %removes states associated with this generator.
            rm_states(1,i+1:end) = rm_states(1,i+1:end) + length(s_rm);
            if ~isempty(qp.Organize.Equalities) && qp.Organize.Equalities(i,1)>0%remove any associated equality constraints
                eq = qp.Organize.Equalities(i,1) + (0:(qp.Organize.Equalities(i,2)-1));
                req_keep(eq) = false;%equality constraints
                qp.Organize.Equalities(i,1:2) = [0,0];
            end
            if qp.Organize.SpinRow(i)>0
                req_keep(qp.Organize.SpinRow(i)) = false;% spinning reserve equality constraint
            end
            if ~isempty(qp.Organize.Inequalities) && qp.Organize.Inequalities(i)>0
                r_keep(qp.Organize.Inequalities(i)) = false;%inequality constraints
            end
        else
            if ~isempty(qp.organize{1,i})%don't remove, but change state #'s
                qp.organize{1,i} = qp.organize{1,i} - rm_states(1,i); %
            end
        end
    end
    %%keep track of states that will be pulled out in sort_solutions_step
    if ~isempty(qp.Organize.HeatVented)
        for i = 1:1:length(qp.Organize.HeatVented)
            qp.Organize.HeatVented(i) = qp.Organize.HeatVented(i) - rm_states(n_g); 
        end
    end
    if ~isempty(qp.Organize.CoolVented)
        for i = 1:1:length(qp.Organize.CoolVented)
            qp.Organize.CoolVented(i) = qp.Organize.CoolVented(i) - rm_states(n_g); 
        end
    end
    if ~isempty(qp.Organize.SpinReserveStates) && any(qp.Organize.SpinReserveStates)
        for i = 1:1:length(qp.Organize.SpinReserveStates(1,:))
            if qp.Organize.SpinReserveStates(i)>0
                qp.Organize.SpinReserveStates(i) = qp.Organize.SpinReserveStates(i) - rm_states(1,min(i,n_g));
            end
        end
    end
%     %keep track of constraint rows, probably not necessary
%     for i = 1:1:n_g
%         if ~isempty(qp.Organize.Equalities) && qp.Organize.Equalities(i,1)>0
%             qp.Organize.Equalities(i,1) = qp.Organize.Equalities(i,1) - nnz(~req_keep(1:qp.Organize.Equalities(i,1)));
%         end
%         if ~isempty(qp.Organize.fluid_loop) && qp.Organize.fluid_loop(i,1)>0
%             qp.Organize.fluid_loop(i,1) = qp.Organize.fluid_loop(i,1) - nnz(~req_keep(1:qp.Organize.fluid_loop.req(i,1)));
%         end
%         if qp.Organize.SpinRow(i)>0
%             qp.Organize.SpinRow(i) = qp.Organize.SpinRow(i) - nnz(~req_keep(1:qp.Organize.SpinRow(i)));
%         end
%         if ~isempty(qp.Organize.Inequalities) && qp.Organize.Inequalities(i)>0
%             qp.Organize.Inequalities(i) = qp.Organize.Inequalities(i) - nnz(~r_keep(1:qp.Organize.Inequalities(i)));
%         end
%     end
    %%eliminate states and constraints to ensure minimum problem realization
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
end%Ends function disable_generators_step