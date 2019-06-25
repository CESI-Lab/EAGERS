function [cost,elec_imbalance,heat_imbalance] = net_cost(gen,var1,date,method)
letters = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';};
if strcmp(method,'Dispatch')
    dispatch = var1;
    input = 0*dispatch;
    elec_imbalance = 0*date;
    heat_imbalance = 0*date;
    for i = 1:1:length(gen)
        skip = false;
        if ~isempty(gen(i).Output)
            cap = gen(i).Output.Capacity*gen(i).Size;
            cap(end) = cap(end)+1e-10; %prevent Nan's during interpolation
        end
        if strcmp(gen(i).Type,'Electric Generator')
            if isfield(gen(i).Output,'Electricity')
                eff = gen(i).Output.Electricity;
            elseif isfield(gen(i).Output,'DirectCurrent')
                eff = gen(i).Output.DirectCurrent;
            end
        elseif strcmp(gen(i).Type,'CHP Generator')
            if isfield(gen(i).Output,'Electricity')
                eff = gen(i).Output.Electricity;
            elseif isfield(gen(i).Output,'DirectCurrent')
                eff = gen(i).Output.DirectCurrent;
            end
            heat = gen(i).Output.Heat;
            CHP = interp1(cap,heat,dispatch(:,i));
            heat2 = dispatch(:,i)./interp1(cap,eff,dispatch(:,i)).*CHP;
            h_gen = -gen(i).QPform.constDemand.H*(dispatch(:,i)>0);
            output = dispatch(:,i);
            n = dispatch(:,i)>0;
            for j = 1:1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end)))
                h_gen = h_gen + min(output,gen(i).QPform.(letters{j}).ub(end))*gen(i).QPform.output.H(j,end);
                output = max(0,output - min(output,gen(i).QPform.(letters{j}).ub(end)));
            end
            heat_imbalance(n) = heat_imbalance(n) + (heat2(n) - h_gen(n));
        elseif strcmp(gen(i).Type,'Chiller') 
            eff = gen(i).Output.Cooling;
            cop = interp1(cap,eff,dispatch(:,i));
            cop(isnan(cop)) = 1;
            skip = true;
            n = dispatch(:,i)>0;
            if strcmp(gen(i).Source,'Electricity')
                elec_use = gen(i).QPform.constDemand.E*(dispatch(:,i)>0);
                output = dispatch(:,i);
                for j = 1:1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end)))
                    elec_use = elec_use - min(output,gen(i).QPform.(letters{j}).ub(end))*gen(i).QPform.output.E(j,end);
                    output = max(0,output - min(output,gen(i).QPform.(letters{j}).ub(end)));
                end
                elec_imbalance(n) = elec_imbalance(n) + (dispatch(n,i)./cop(n) - elec_use(n));
            elseif strcmp(gen(i).Source,'Heat')
                heat_use = gen(i).QPform.constDemand.H*(dispatch(:,i)>0);
                output = dispatch(:,i);
                for j = 1:1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end)))
                    heat_use = heat_use - min(output,gen(i).QPform.(letters{j}).ub(end))*gen(i).QPform.output.H(j,end);
                    output = max(0,output - min(output,gen(i).QPform.(letters{j}).ub(end)));
                end
                heat_imbalance(n) = heat_imbalance(n) + (dispatch(n,i)./cop(n) - heat_use(n));
            end
        elseif strcmp(gen(i).Type,'Heater')
            eff = gen(i).Output.Heat;  
        else
            skip = true;
        end
        if ~skip %dont add the cost of a chiller if you ran E and C simultaneously, or you will double count the chiller demand
            input(:,i) = dispatch(:,i)./interp1(cap,eff,dispatch(:,i));
            if nnz(dispatch(:,i))<length(date) %if you have steps where the dispatch is zero, prevent Nan
                OffSteps = (dispatch(:,i)==0);
                input(OffSteps,i) = 0;
            end
        elseif strcmp(gen(i).Type,'Utility')
            input(:,i) = dispatch(:,i);
        end
    end
    %% ___ %%% alternate: uses the cost in the optimization (FitB)
%     dt = Timestamp(2:end)-Timestamp(1:end-1);
%     Cost = zeros(length(Timestamp)-1,1);
%     StorPower =[];
%     for i = 1:1:length(Gen)
%         if nnz(Dispatch(:,i))>0
%             if ~isempty(Gen(i).QPform.states)
%                 states = Gen(i).QPform.states(:,end);
%             end
%             if isfield(Gen(i).QPform,'Stor')
%                 StorPower(:, end+1) = Dispatch(1:end-1,i) - Dispatch(2:end,i); %no cost
%             elseif isfield(Gen(i).QPform,'constCost') %all of these cost terms need to be scaled later on
%                 I = Gen(i).QPform.(states{1}).ub(end);
%                 b1 = Gen(i).QPform.(states{1}).f(end);
%                 b2 = Gen(i).QPform.(states{2}).f(end);
%                 b3 = Gen(i).QPform.(states{2}).H(end);
%                 b4 = Gen(i).QPform.constCost;
%                 Cost = Cost + (min(Dispatch(2:end,i),I)*b1 + max(Dispatch(2:end,i)-I,0)*b2 + max(Dispatch(2:end,i)-I,0).^2*b3+b4).*dt.*scaleCost(:,i);
%             elseif~isempty(states) %single state generators (linear cost term) & utilities
%                 b1 = Gen(i).QPform.(states{1}).f(end);
%                 Cost = Cost + Dispatch(2:end,i)*b1.*scaleCost(:,i).*dt;
%                 %need to handle case of sell-back
%             elseif ~isempty(strcmp(Gen(i).Source,'Renewable'))
%                 %renewable
%             end
%         end
%     end  
elseif strcmp(method,'Input')
    input = var1;
    dispatch = 0*input;
    dispatch(input>1) = 1;
end
n_s = nnz(date);
n_g = length(gen);
dt = (date(2:n_s) - date(1:n_s-1))*24;
scale_cost = update_cost(date(1:n_s-1),gen).*(dt*ones(1,n_g));%% All costs were assumed to be 1 when building matrices
stor =[];
startup_cost = zeros(n_s-1,n_g);
for i = 1:1:n_g
    if isfield(gen(i).VariableStruct, 'StartCost')
        n_startups = (dispatch(2:n_s,i)>0).*(dispatch(1:n_s-1,i)==0);
        startup_cost(:,i) = n_startups*gen(i).VariableStruct.StartCost;
    end
    if isfield(gen(i).QPform,'Stor')
       stor(end+1) = i; 
    end
end
scale_cost(:,stor)=0; %remove cost of storage
cost = sum(input(2:n_s,1:n_g).*scale_cost + startup_cost,2);