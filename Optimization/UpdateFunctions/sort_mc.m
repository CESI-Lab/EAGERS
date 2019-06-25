function sort_mc = sort_mc(margin_cost,out,t_max,off_gen,dt)
%Converts MarginCost array to 5 columns, sorting by the lowest marginal 
%generation block cost before time t 
%1 = cumulative marginal energy, 
%2 = cumulative marginal cost, used for interpolating to see which blocks are used
%3 = Generator that increases output, 
%4 = time point at which generator increases output, 
%5 = amount that generator increases output (kW)
[n_g,n_s,n] = size(margin_cost.(out).Capacity.SpinReserve);
mc = zeros(n_g*n_s*n+1,4);%capacity, cost per kW, generator, time
k = 1;
for i = 1:1:n_g
    if isempty(off_gen) || i~=off_gen
        for t = 1:1:t_max
            if margin_cost.(out).Capacity.SpinReserve(i,t,1)>0
                for j = 1:1:n
                    mc(k+j,1) = margin_cost.(out).Capacity.SpinReserve(i,t,j);
                    mc(k+j,2) = margin_cost.(out).Cost.SpinReserve(i,t,j);
                end
                mc(k+1:k+n,3) = i;
                mc(k+1:k+n,4) = t;
                mc(k+1:k+n,5) = dt(t);
                k = k+n;
            end
        end
    end
end
if k>1
    mc = mc(2:k,:);
    sort_mc = zeros(k-1,5);
    [~,I] = sort(mc(:,2));
    sort_mc(:,3) = mc(I,3);
    sort_mc(:,4) = mc(I,4);
    sort_mc(:,5) = mc(I,1);
    sort_mc(1,1) = sort_mc(1,5)*mc(I(1),5);% (kW*dt) cumulative energy instead of power
    sort_mc(1,2) = sort_mc(1,5)*mc(I(1),2)*mc(I(1),5);% (kW*$/kW *dt) cumulative cost of the net kJ
    for r = 2:1:k-1
        sort_mc(r,1) = sort_mc(r-1,1) + sort_mc(r,5)*mc(I(r),5);
        sort_mc(r,2) = sort_mc(r-1,2) + sort_mc(r,5)*mc(I(r),2)*mc(I(r),5);%cumulative cost = previous total + size*cost per kW * time
    end
else
    sort_mc = [];
end
end%ends function sort_mc