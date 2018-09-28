function solution = fix_infeasible(qp)
[m,n] = size(qp.organize);
n_s = max(1,m-1);
s_add = 0;
f_name = fieldnames(qp.Organize.Balance);
for f = 1:1:length(f_name)
    s_add = s_add + n_s*length(qp.Organize.Balance.(f_name{f})(1,:));
end

[r,~] = size(qp.A);
[req,x_l] = size(qp.Aeq);
qp.H = diag([diag(qp.H);zeros(s_add,1);]);
qp.f = [qp.f;1e3*max(qp.f)*ones(s_add,1);];
qp.Aeq = [qp.Aeq,zeros(req,s_add);];
if r>0 
    qp.A = [qp.A,zeros(r,s_add)];
end
qp.lb = [qp.lb;0*ones(s_add,1);];
qp.ub = [qp.ub;1e10*ones(s_add,1);];  

k = 0;
for f = 1:1:length(f_name)
    for t = 1:1:n_s
        for j = 1:1:length(qp.Organize.Balance.(f_name{f})(1,:))
            k = k+1;
            qp.Aeq(qp.Organize.Balance.(f_name{f})(t,j),x_l+k) = 1;
        end
    end
end

[x,feasible] = call_solver(qp);
if feasible ~=1
    solution = x;
    disp('some additional contraints making problem infeasible, tried adding infinite sources');
else
    solution = sort_solution(x,qp);
    k = 0;
    for f = 1:1:length(f_name)
        j = length(qp.Organize.Balance.(f_name{f})(1,:));
        n_makeup = nonzeros((1:(n_s*j))'.*(x(x_l+k+1:x_l+k+n_s*j)>1e-3));
        k = k+n_s*j;
        if ~isempty(n_makeup)
            if length(n_makeup)<=4
                for i = 1:1:length(n_makeup)
                    disp(strcat('Short production capacity of ',f_name{f}, ' at node ',num2str(rem(n_makeup(i),j)),' at time intevals ',num2str(ceil(n_makeup(i)/j))))
                end
            else
                disp(strcat('Short production capacity of ',f_name{f}, ' at ',num2str(length(n_makeup)),' time intevals or nodes'))
            end
        end
    end
end
end%ends funtion fix_infeasible