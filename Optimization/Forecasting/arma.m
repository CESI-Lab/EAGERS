function forecast = arma(date,prev)
%% Need to revamp
forecast.Timestamp = date;
alpha = .85;
beta = 0.9; % Take this out after training is added!
outs = fieldnames(prev.Demand);
for S = 1:1:length(outs)
    if ~isfield(prev,'prevErr') || ~isfield(prev.prevErr,outs{S}) || isnan(prev.prevFcast.(outs{S}))
        prev.prevErr.(outs{S}) = 0;
        prev.prevFcast.(outs{S}) = prev.Demand.(outs{S});
    else
        prev.prevErr.(outs{S}) = prev.prevFcast.(outs{S}) - prev.Demand.(outs{S})(end) + prev.prevErr.(outs{S});
    end
end
d_vec = prev.Timestamp;
n_days = ceil(date(end)-prev.Timestamp(end));
w = [linspace(alpha,0,nnz(date<(prev.Timestamp(end)+beta)))';zeros(nnz(date>=(prev.Timestamp(end)+beta)),1)];
n = length(d_vec)-1;
for d = 1:1:n_days
    d_vec(end+1:end+n) = prev.Timestamp(2:end) + d;
end
for S = 1:1:length(outs)
    data = prev.Demand.(outs{S});
    for d = 1:1:n_days
        data(end+1:end+n) = prev.Demand.(outs{S})(2:end);
    end
    forecast.Demand.(outs{S}) = interp1(d_vec,data,date);
    forecast.Demand.(outs{S})(1:length(w)) = forecast.Demand.(outs{S})(1:length(w)) + w*prev.prevErr.(outs{S}); 
end
end%Ends function arma