function [SolutionGurobi,Solution_mcQP,Solution_cQP, fitADispatch, tsim] = GurobiTest(Date)
global Plant TestData
load_test_data
Date = Date+[0;build_time_vector(Plant.optimoptions)/24];

if ~isfield(Plant,'Building')
    Plant.Building = [];
end
if ~isfield(Plant,'fluid_loop') 
    Plant.fluid_loop = [];
end
if isfield(Plant,'Data') 
    data = Plant.Data;
else
    data = [];
end
if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
    [TestData,Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.subNet,Plant.OpMatA,Plant.OpMatB,Plant.OneStep,Plant.Online] = initialize_optimization(Plant.Generator,Plant.Building,Plant.fluid_loop,Plant.Network,Plant.optimoptions,data,TestData);
end
Buildings = Plant.Building;
fluid_loop = Plant.fluid_loop;

tsim = zeros(1,5);

TestData.RealTimeData = interpolate_data(TestData,Plant.optimoptions.Resolution*3600,0.00);%create test data at correct frequency
Plant.Generator = automatic_ic(Plant.Generator,Buildings,fluid_loop,Plant.subNet,Date,Plant.OneStep,Plant.optimoptions);% set the initial conditions
[Forecast,Plant.Generator,Buildings] = update_forecast(Plant.Generator,Buildings,fluid_loop,Plant.subNet,Plant.optimoptions,Date(2:end));
dt = (Date(2:end) - Date(1:end-1))*24;
nG = length(Plant.Generator);
nS = length(Date)-1;
scaleCost = update_cost(Date(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
PredictDispatch = zeros(length(Date),nG);
for i = 1:1:nG
    PredictDispatch(:,i) = Plant.Generator(i).CurrentState;
end
marginCost = update_mc(Plant.Generator,PredictDispatch,scaleCost,dt,[]);
Plant.optimoptions.MixedInteger = true;
QP = update_matrices(Plant.Generator,Buildings,fluid_loop,Plant.subNet,Plant.optimoptions,Plant.OpMatA,Date,scaleCost,marginCost,Forecast,[]);
%% Step 1 Determine initial dispatch
Locked = true(nS+1,nG);
for i = 1:1:nG
    if ~Plant.Generator(i).Enabled
        Locked(:,i) = 0;
    end
end
tic
QP = disable_generators(QP,Locked,[]);%Disable generators here
[x,Feasible] = call_solver(QP);
if Feasible == 1
    Solution = sort_solution(x,QP);
end
tsim(1,2) = toc;
if ~(Feasible==1)%% hopefully not here
    disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
    Solution.Dispatch = PredictDispatch;
end
fitADispatch = Solution.Dispatch;

%% Step 2:  dispatch step by step
tic
OptimalState = dispatch_step(Plant.Generator,Buildings,fluid_loop,Plant.subNet,Plant.optimoptions,Plant.OneStep,Date(1),Forecast,scaleCost,dt,Solution.Dispatch);
clear mex
tsim(1,2) = toc;
%% Start with optimal dispatch, and check if feasible
tic
marginCost = update_mc(Plant.Generator,OptimalState,scaleCost,dt,[]);
QP_0 = update_matrices(Plant.Generator,Buildings,fluid_loop,Plant.subNet,Plant.optimoptions,Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
Locked = verify_ramping(Plant.Generator,Plant.subNet,QP_0,Locked,OptimalState,dt);
QP = disable_generators(QP_0,Locked,[]);%Disable generators here
[x,Feasible] = call_solver(QP);%this is the dispatch with fit B
if Feasible == 1
    Solution_mcQP = sort_solution(x,QP);
else
    Solution_mcQP = [];
end
tsim(1,4) = toc;


%% Run Gurobi
PredictDispatch = zeros(length(Date),nG);
for i = 1:1:nG
    PredictDispatch(:,i) = Plant.Generator(i).CurrentState;
end
marginCost = update_mc(Plant.Generator,PredictDispatch,scaleCost,dt,[]);
QP = update_matrices(Plant.Generator,Buildings,fluid_loop,Plant.subNet,Plant.optimoptions,Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]);
QP.constCost = dt*QP.constCost.*scaleCost;
nG = length(Plant.Generator);

%add in end conditions so that storage is the same at the final timestep
for i = 1:1:nG
    if isfield(Plant.Generator(i).QPform,'Stor')
        [m,n] = size(QP.Aeq);
        QP.Aeq = [QP.Aeq;zeros(1,n);];
        QP.Aeq(m+1,QP.organize{end,i}) = 1;
        QP.beq = [QP.beq;Solution_mcQP.Dispatch(end,i);];
    end
end
tic
x = gurobi_opt(QP);
if ~any(isnan(x))
    SolutionGurobi = sort_solution(x,QP);
else 
    SolutionGurobi = [];
end
tsim(1,1) = toc;


%% run cQP 
tic
marginCost = update_mc(Plant.Generator,Solution1.Dispatch,scaleCost,dt,[]);
QP_0 = update_matrices(Plant.Generator,Buildings,fluid_loop,Plant.subNet,Plant.optimoptions,Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
[Solution_cQP,LBrelax]= cqp_method(Plant.Generator,QP_0,Solution.Dispatch,Date);
tsim(1,5) = toc;


%% compare costs
if ~isempty(Solution_cQP)
    [Cost,Eimbalance,Himbalance] = net_cost(Plant.Generator,Solution_cQP.Dispatch,Date,'Dispatch');
    Solution_cQP.Cost = Cost;
    Solution_cQP.Eimbalance = Eimbalance;
    Solution_cQP.Himbalance = Himbalance;
    Solution_cQP.LBrelax = LBrelax;
end

if ~isempty(Solution_mcQP)
    [Cost,Eimbalance,Himbalance] = net_cost(Plant.Generator,Solution_mcQP.Dispatch,Date,'Dispatch');
    Solution_mcQP.Cost = Cost;
    Solution_mcQP.Eimbalance = Eimbalance;
    Solution_mcQP.Himbalance = Himbalance;
end

if ~isempty(SolutionGurobi)
    [Cost,Eimbalance,Himbalance] = net_cost(Plant.Generator,SolutionGurobi.Dispatch,Date,'Dispatch');
    SolutionGurobi.Cost = Cost;
    SolutionGurobi.Eimbalance = Eimbalance;
    SolutionGurobi.Himbalance = Himbalance;
end

