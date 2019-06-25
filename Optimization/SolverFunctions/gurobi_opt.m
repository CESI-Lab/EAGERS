function [x,Feasible] = gurobi_opt(QP)
global Plant
if ~isfield(QP,'Organize') %calculating generator cost fits
    n = length(QP.f);
    if ~isempty (QP.lb)
        QP.A = [QP.A;-eye(n);];
        QP.b = [QP.b;QP.lb;];
    end
    if ~isempty (QP.ub)
        QP.ub(isinf(QP.ub)) = 1e3*max(QP.ub(~isinf(QP.ub)));
        QP.A = [QP.A;eye(n);];
        QP.b = [QP.b;QP.ub;];
    end
    cvx_begin quiet
        variable x(n) 
        minimize (0.5*x'*QP.H*x+QP.f'*x)
        subject to
            QP.Aeq*x == QP.beq;
            QP.A*x <= QP.b;
    cvx_end
else %solving mixed-integer generator dispatch
    nG = length(Plant.Generator);
    nS = length(QP.organize(:,1))-1;
    ic = nnz(QP.Organize.IC);
    totalStates = nS*QP.Organize.t1States+ic;
    disp = find(QP.Organize.Dispatchable>0);
    nDG = length(disp);
    boolVar = nS*nDG;
    D = zeros(1,boolVar);
    for k = 1:1:nDG
        D(1,k:nDG:boolVar) = QP.constCost(:,disp(k))';
    end
    H = full(QP.H);
    f = QP.f;
    F = fieldnames(QP.constDemand);
    %organize constant demands
    constDem = zeros(length(QP.beq),boolVar);
    for i = 1:1:length(F)
        if strcmp(F{i},'E')
            net = 'Electrical';
        elseif strcmp(F{i},'H')
            net = 'DistrictHeat';
        elseif strcmp(F{i},'C')
            net = 'DistrictCool';
        end
        for n = 1:1:length(Plant.subNet.(net).nodes) %run through all the nodes in this network
            equip = Plant.subNet.(net).Equipment{n}; %equipment at this node
            req = QP.constDemand.(F{i}).req((n-1)*nS+1);
            for j = 1:1:length(equip)
                if any(equip(j)==disp)
                    k = equip(j);
                    if isfield(Plant.Generator(k).QPform,'constDemand')  && isfield(Plant.Generator(k).QPform.constDemand,F{i})
                        boolIndex = min(nonzeros((1:1:nDG).*(equip(j)==disp)));
                        for t = 1:1:nS
                            constDem(req+(t-1)*QP.Organize.t1Balances,boolIndex+(t-1)*nDG) = Plant.Generator(k).QPform.constDemand.(F{i});
                        end
                    end
                end
            end
        end
    end
    %organize which states are dispatchable
    UB = zeros(totalStates,boolVar+1);
    LB = zeros(totalStates,boolVar+1);
    UB(:,end) = QP.ub;
    LB(:,end) = QP.lb;
    k = 0;
    for i = 1:1:nG
        states = QP.Organize.States{i};
        if QP.Organize.Dispatchable(i)
            k = k + 1;
            for t = 1:1:nS %move upper & lower constrain to multiply by boolean
                UB(states+(t-1)*QP.Organize.t1States,k+(t-1)*nDG) = QP.ub(states+(t-1)*QP.Organize.t1States);
                LB(states+(t-1)*QP.Organize.t1States,k+(t-1)*nDG) = QP.lb(states+(t-1)*QP.Organize.t1States);
                UB(states+(t-1)*QP.Organize.t1States,end) = 0;
                LB(states+(t-1)*QP.Organize.t1States,end) = 0;
            end
        end
    end

    cvx_begin %quiet
    cvx_solver_settings('timelimit', 3600)
    % cvx_solver_settings('presolve', 0)
    % cvx_solver_settings( 'dumpfile', 'test' ) 
    % cvx_solver_settings( 'NumericFocus', '3' ) 

        variable x(totalStates)
        variable q(boolVar) binary

        C = x'*H*x + f'*x + D*q;

        minimize ( C )

        subject to
            QP.A*x <= QP.b;%all inequality constraints
            QP.Aeq*x == QP.beq + constDem*q;
            LB*[q;1] <= x <= UB*[q;1];
    cvx_end
end
if strcmp(cvx_status,'Solved')
    Feasible = 1;
else
    Feasible = -1;
end
end%Ends function gurobi_opt