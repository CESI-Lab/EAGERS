function z = typical_day(h,timestamp,data)
steps = round(1/(timestamp(2)-timestamp(1)));%points in 1 day of data
resolution =24/steps;
[~,m,~,hour,minutes] = datevec(timestamp(2:end-1));%avoid checking month 1st and last point if starting or ending at 00
m = [m(1); m; m(end)];%make m the same length as data
hour = [max(0,floor(hour(1)-resolution)); hour; min(24,hour(end)+floor((minutes(end)+60*resolution)/60))];%make hour the same length as data
months = sort(unique(m));
%% fit temperature data
z = zeros(12,24);
for i = 1:1:length(months)
    if ~isempty(h)
        waitbar(i/length(months),h,'Calculating Monthly Profiles');
    end
    for j = 1:1:24
        X = (m==i) & (hour == j-1);
        if any(data(X)~=0)
            z(months(i),j) = sum(data(X))./nnz(X);
        end
    end
end
if length(months)<12
    if nnz(months==1)==1 && nnz(months==12)==1
        j =1;
        for i= 2:1:11
            if nnz(months==i)==1
                j=j+1;
            else z(i,:) = z(j,:);
            end
        end
    else
        for i= 1:1:min(months)-1
            z(i,:) = z(min(months),:);
        end
        for i= max(months)+1:1:12
            z(i,:) = z(max(months),:);
        end
    end
end
end%Ends function typical_day