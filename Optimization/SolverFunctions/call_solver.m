function [x,feasible] = call_solver(qp)
if strcmp(qp.solver,'ANN')
    qp.solver = 'quadprog';
end
if strcmp(qp.solver,'linprog') && any(any(qp.H - diag(diag(qp.H)))) %not a seperable QP problem
    qp.solver = 'quadprog';
end
if strcmp(qp.solver,'quadprog') && ~license('test','Optimization_Toolbox')
    qp.solver = 'qpip';
end
if isempty(qp.f)
    x = [];
    feasible = 1;
else
    switch qp.solver
        case 'quadprog'
        %use matlabs linprog or quadprog
        if ~any(qp.H)
            options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none'); %,'Display','iter-detailed');% ,'TolFun',1e-10,'TolX',1e-10
            [x,~,feasible] = linprog(qp.f,qp.A,qp.b,qp.Aeq,qp.beq,qp.lb,qp.ub,[],options); 
        else
            options = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');%,'TolFun',1e-10,'TolX',1e-10
            [x,~,feasible] = quadprog(qp.H,qp.f,qp.A,qp.b,qp.Aeq,qp.beq,qp.lb,qp.ub,[],options);
        end
        case 'Gurobi'
            [x,feasible] = gurobi_opt(qp);
        case 'qpip'
            [x,feasible] = qpip(qp);
    end
end
end%Ends call_solver