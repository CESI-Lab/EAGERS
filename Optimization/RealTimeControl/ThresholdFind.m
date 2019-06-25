function [Feasible,Cost]= ThresholdFind(Output,Demand,DemandH,QPall,Organize,constCost,EnabledA,EnabledB)
options = optimset('MaxIter',50,'Display','none');
options2 = optimset('Algorithm','interior-point-convex','MaxIter',50,'Display','none');
%% ----- %%%% %find demand threshold where changeover should occur
%% Switch from EnabledA to EnabledB at time steps 1:n
%% find infeasible cases, and best cost case as forecasted demand rises/falls
[n,r] = size(Demand);
Cost = zeros(n,n);
Feasible = zeros(n,n);
Enable = zeros(n,length(EnabledA));
for k = 1:1:r % change energy demand within range
    %update beq in QPall
    mat = Organize.(Output).Demand{1};
    index = Organize.(Output).Demand{2};
    QPall.(mat)(index) = Demand(:,k);
    if ~isempty(DemandH) %put heat in with electricity
        mat = Organize.H.Demand{1};
        index = Organize.H.Demand{2};
        QPall.(mat)(index) = DemandH(:,k);
    end
    for t = 1:1:n %change when swichover occurs
        for i = 1:1:t
            Enable(i,:) = EnabledA;
        end
        for i = t+1:1:n
            Enable(i,:) = EnabledB;
        end
        QP = disableGenerators(QPall,Organize,Enable,[]);%Disable generators here
        if nnz(QP.H)==0
            [~,cost,Feas] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
        else [~,cost,Feas] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2); 
        end
        if Feas~=1
            Feasible(k,t) = NaN;
            Cost(k,t) = inf;
        else Cost(k,t) = cost + sum(constCost.*sum(Enable)); 
            Feasible(k,t) = 1;
        end
    end
end
