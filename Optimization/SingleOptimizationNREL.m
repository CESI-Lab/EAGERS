function solution = SingleOptimizationNREL(date1)
global Plant TestData
n_g = length(Plant.Generator);%skip initialization
ic = zeros(1,n_g);
for i=1:1:n_g
    if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
        ic(i) = 0.5*Plant.Generator(i).Size*(Plant.Generator(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
    end
end
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    Buildings = Plant.Building;
else
    Buildings = [];
end

if isfield(Plant,'Data') 
    data = Plant.Data;
else
    data = [];
end
TestData = update_test_data(TestData,data,Plant.Generator,Plant.optimoptions);
date = date1+[0;build_time_vector(Plant.optimoptions)/24];%linspace(Date,DateEnd)';would need to re-do optimization matrices for this time vector
freq = 1; %period of repetition (1 = 1 day)
res = Plant.optimoptions.Resolution/24;
n_o = round(freq/res)+1;
prev_data = get_data(TestData,linspace((date1 - res - freq),date1-res,n_o)',[],[]);
now_data = get_data(TestData,date1,[],[]);
future_data = get_data(TestData,date(2:end),[],[]);

[Forecast,Plant.Generator,~] = update_forecast(Plant.Generator,Buildings,[],[],Plant.optimoptions,date(2:end),TestData.HistProf,prev_data,now_data,future_data);
scaleCost = update_cost(date(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
solution = solver_nrel(Plant.Generator,Plant.optimoptions,ic(1:n_g),22,Forecast,scaleCost);
solution.Timestamp = date;
Plant.Building.Tzone = solution.Buildings.Temperature(1);
end %Ends function SingleOptimizationNREL