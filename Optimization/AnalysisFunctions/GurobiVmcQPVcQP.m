function GurobiVmcQPVcQP(Date0, Date_end)
%% GurobiVmcQPVcQP runs and saves the results from GurobiTest, which 
% compares the dispatches from Gurobi, mcQP, and cQP methods. the solutions
% are saved in the folder GUI\Optimization\Results
% INPUTS:
%   Date0       -   The date vector for the start time for the run in the
%                   format [YYYY, M, D, H, M, S] don't include zeros
%   Date_end    -   The date vector for the last dispatch for the run in
%                   the format same as above
% 
% load Plant file from directory before running
global Model_dir Plant


FirstDay = datenum(Date0);
LastDay = datenum(Date_end);
dt = Plant.optimoptions.Resolution;%hourly timesteps
Dates = FirstDay:dt/24:LastDay;
nD = length(Dates);%number of days
nS = Plant.optimoptions.Horizon/dt;%number of steps
nG = length(Plant.Generator);%number of components

%preallocate space
allSols = [];
allSols.Timers = zeros(nD,5);%gurobi, fit A, mcQP step2, mcQP step3, cQP step2
allSols.Timestamp = Dates;
allSols.Cost = zeros(nS,3,nD);
allSols.Eimbalance = zeros(nS+1,3,nD);
allSols.Himbalance = zeros(nS+1,3,nD);
allSols.Dispatch = [];
allSols.Dispatch.Gurobi = zeros(nS+1,nG,nD);
allSols.Dispatch.mcQP = zeros(nS+1,nG,nD);
allSols.Dispatch.cQP = zeros(nS+1,nG,nD);
allSols.Dispatch.fitA = zeros(nS+1,nG,nD);
allSols.cQPLBrelax = zeros(nD);

for Si = 1:1:nD
    Date = Dates(Si);
    [GurobiSol, mcQPSol, cQPSol, fitAsol, Timers] = GurobiTest(Date);
    
    allSols.Timers(Si,:) = Timers;
    
    if ~isempty(GurobiSol)%only try to save if it was feasible solution
        allSols.Dispatch.Gurobi(:,:,Si) = GurobiSol.Dispatch;
        allSols.Cost(:,1,Si) = GurobiSol.Cost;
        allSols.Eimbalance(:,1,Si) = GurobiSol.Eimbalance;
        allSols.Himbalance(:,1,Si) = GurobiSol.Himbalance;
    end
    
    if ~isempty(mcQPSol)
        allSols.Dispatch.mcQP(:,:,Si) = mcQPSol.Dispatch;
        allSols.Cost(:,2,Si) = mcQPSol.Cost;
        allSols.Eimbalance(:,2,Si) = mcQPSol.Eimbalance;
        allSols.Himbalance(:,2,Si) = mcQPSol.Himbalance;
    end
    
    if ~isempty(cQPSol)
        allSols.Dispatch.cQP(:,:,Si) = cQPSol.Dispatch;
        allSols.Cost(:,3,Si) = cQPSol.Cost; 
        allSols.Eimbalance(:,3,Si) = cQPSol.Eimbalance;
        allSols.Himbalance(:,3,Si) = cQPSol.Himbalance;
        allSols.cQPLBrelax(Si) = cQPSol.LBrelax;
    end
    
    allSols.Disptach.fitA(Si,:,:) = fitAsol;
   
    %save periodically
    if rem(Si,1)==0
        save(fullfile(Model_dir,'GUI\Optimization\Results',strcat(Plant.Name,'allSols',num2str(FirstDay),'.mat')),'allSols')
    end
end
    
%save at the end
save(fullfile(Model_dir,'GUI\Optimization\Results',strcat(Plant.Name,'allSols',num2str(FirstDay),'.mat')),'allSols')
    
    
    
    
    
