function test_data = update_test_data(test_data,data,weather,options)
if ~isfield(test_data,'Demand') && isfield(data,'Demand')
    test_data.Demand = data.Demand;
end
if ~isfield(test_data,'Weather') && ~isempty(weather)
    test_data.Weather = weather;
elseif ~isfield(test_data,'Weather') && isfield(data,'Weather')
    test_data.Weather = data.Weather;
else
    %load weather
end
if isfield(data,'Hydro') && ~isfield(test_data,'Hydro')
    test_data.Hydro = data.Hydro;
end
if ~isfield(test_data,'HistProf')
    test_data.HistProf = [];
end
if isempty(test_data.HistProf) &&  isfield(data,'HistProf') && ~isempty(data.HistProf)
    test_data.HistProf = data.HistProf;
end
%create typical day fits if necessary
if ~strcmp(options.forecast,'Perfect') && isempty(test_data.HistProf)
    s = fieldnames(test_data.Weather);
    for j = 1:1:length(s)
        if isnumeric(test_data.Weather.(s{j}))
            test_data.HistProf.(s{j}) = typical_day([],test_data.Timestamp,test_data.Weather.(s{j}));
        end
    end
    if isfield(test_data,'Hydro')
        nodes = length(test_data.Hydro.SourceSink(1,:));
        for n = 1:1:nodes
            test_data.HistProf.Hydro.SourceSink(n) = {typical_day([],test_data.Hydro.Timestamp,test_data.Hydro.SourceSink(:,n))};
        end
    end
end
%create surface fits if necessary
if strcmp(options.forecast,'Surface') && isfield(test_data,'Demand') 
    s = fieldnames(test_data.Demand);
    if ~isfield(test_data.HistProf,'Demand') || ~all(ismember(s,fieldnames(test_data.HistProf.Demand)))
        test_data.HistProf.Demand = calculate_fits(test_data.Demand,test_data.Timestamp,test_data.Weather);%% calculate surface fits used in forecasting
    end
end
end%ends function update_test_data