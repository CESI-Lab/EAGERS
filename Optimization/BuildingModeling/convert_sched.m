function newsched = convert_sched(sched,ramp)
[m,n] = size(sched);
if m==2
    newsched = sched; %this is constant all day, already made 0 hour and 24 hour
else
    newsched = zeros(2*(m-1),n);
    sched(1,2) = sched(2,2);
    newsched(1,:) = sched(1,:);%hour 0
    newsched(end,:) = sched(m,:);%hour 24
    if ramp<1e-3
        for i = 2:1:m-1 
            newsched(2*i-2,1) = sched(i);
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+ramp;
            newsched(2*i-1,2) = sched(i+1,2);
        end
    else 
        for i = 2:1:m-1 %add points in the middle so it can be properly interpolated
            t_bef = sched(i) - sched(i-1);
            t_aft = sched(i+1)-sched(i);
            newsched(2*i-2,1) = sched(i)-min(ramp,t_bef/2);
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+min(ramp,t_aft/2);
            newsched(2*i-1,2) = sched(i+1,2);
        end
    end
end
end%Ends function convert_sched