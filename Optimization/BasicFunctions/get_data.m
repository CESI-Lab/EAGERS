function select_data = get_data(data,date,timestamp,variability)
%load current, past data or future data
if isnumeric(data)
    n_s = length(date);
    index = nnz(timestamp<date(1)+1e-6)+ (0:n_s-1)'; %+1e-6 to avoid rounding problems
    need_interp = true;
    if n_s==1
        need_interp = false;
    else
        if length(timestamp)>=index(end)
            r = (timestamp(index(2:end))-timestamp(index((1:end-1))))./(date(2:end)-date(1:end-1));
            if all(abs(r-1)<1e-7)
                need_interp = false;
            end
        end
    end

    if need_interp && ~isempty(variability)
        noise = create_noise(n_s,variability);%% create noise vector
    else
        noise = ones(n_s,1);
    end
    if ~need_interp%no need to interpolate
        select_data = data(index,:);
    else
        select_data = zeros(n_s,length(data(1,:)));
        x1 = index(1);
        for t = 1:1:n_s
            k = 0;
            r0 = (timestamp(x1+1) - date(t))/(timestamp(x1+1) - timestamp(x1));
            while r0<0 %summing multiple data points
                k = k + 1+r0;
                select_data(t,:) = select_data(t,:) + (1+r0)*data(x1,:);
                x1 = x1+1;
                r0 = (timestamp(x1+1) - date(t))/(timestamp(x1+1) - timestamp(x1));
            end
            if k>0%summing multiple points                
                select_data(t,:) = (select_data(t,:) + (1-r0)*data(x1,:))/(k+(1-r0));
            else %interpolating
                select_data(t,:) =  (r0*data(x1,:) + (1-r0)*data(x1+1,:))*noise(t);
            end
        end
    end
else
    select_data.Timestamp = date;
    if date(1)<data.Timestamp(1)
        date = date + ceil(data.Timestamp(1)-date(1));%add a whole # of days
    end
    x1 = nnz(data.Timestamp<date(1)+1e-6);
    x2 = min(length(data.Timestamp),nnz(data.Timestamp<date(end)+1e-6)+1);
    f = fieldnames(data);
    f = f(~strcmp('Timestamp',f));
    f = f(~strcmp('HistProf',f));
    for j = 1:1:length(f)
        if isstruct(data.(f{j}))
            s = fieldnames(data.(f{j}));
            for i = 1:1:length(s)
                if isnumeric(data.(f{j}).(s{i})) && ~isempty(data.(f{j}).(s{i}))
                    select_data.(f{j}).(s{i}) = get_data(data.(f{j}).(s{i})(x1:x2,:),date,data.Timestamp(x1:x2),variability);
                end
            end
        else
            if isnumeric(data.(f{j})) && ~isempty(data.(f{j}))
                select_data.(f{j}) = get_data(data.(f{j})(x1:x2,:),date,data.Timestamp(x1:x2),variability);
            end
        end
    end 
end
end%Ends function get_data

function noise = create_noise(num_steps,variability)
z = randn(num_steps,1);
noise = zeros(num_steps,1);% scaled noise to a signal with average magnitude of 1
noise(1) = 0;
noise(2) = (rand(1,1)-.5)*variability; %the value .04 makes the final noise signal peaks = variability with this strategy.
b= sign(noise(2)-noise(1));
c = noise(2);
for n = 3:length(noise)
    %if the noise is increasing the probability is such that the noise should continue to increase, but as the noise approaches the peak noise magnitude the probability of switching increases.
    if (c>0 && noise(n-1)>0) || (c<0 && noise(n-1)<0)
        a = 2*(variability - abs(noise(n-1)))/variability; %constant 2 = 97.7% probability load changes in same direction, 1.5 = 93.3%, 1 = 84.4%, .5 = 69.1%
    else
        a = 2;
    end
    if abs(z(n))>.68 %only changes value 50% of the time
         c = b*(z(n)+a); % c is positive if abs(Noise)is increasing, negative if decreasing
         b = sign(c);
         noise(n) = noise(n-1)+c*variability; 
    else
        noise(n) = noise(n-1);
    end
end
end%ends function create_noise