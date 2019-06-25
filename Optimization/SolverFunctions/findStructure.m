function [AccTrain, AccValid, percOn, Tree, TrimmedT] = findStructure()
%this function is used to find the best structure for a neural network to
%be trained by the given information
%the training and validation sets are contained in an optimal dispatch file
%outputs: AccTrain training accuracy rows represent kernel degrees, columns
%represent components, AccValid validation accuracy, Tree is the untrimmed
%decision tree, TrimmedT is the pruned decision tree
%the trees are used to determine which inputs are valuable
%the training and testing accuracy are used to determine the degree needed
%for a good relationship which maps to the number of layers for a neural
%network
global Plant
    ndisps = floor(length(Plant.Predicted.GenDisp(1,1,:)));%*0.001);%number of dispatch samples for training
    nG = length(Plant.Generator);%number of components
    Out = fieldnames(Plant.RunData.Demand);%outputs from the microgrid, eg Electric, Cooling, Heating
    Time = [1:24]';%vector of timesteps
    dt = repmat((Time-[0;Time(1:end-1)]),ndisps,1);
    nS = length(Time); %number of timesteps per horizon
    OnOff = true(ndisps*24,nG);
    GenDisp = zeros(ndisps*25,nG);
    scaleCost = zeros(ndisps*24,nG);
    RenewE = zeros(ndisps*25,1);
    FirstDisp = zeros(ndisps*25,nG);
    %find all dispatchable components
    gens = [];%dispatchable generators
    stor = [];%energy storage units
    renews = [];%renewable sources
    UB = zeros(1,nG);
    LB = zeros(1,nG);
    for i = 1:1:nG
        %find upper and lower bounds so that you can normalize inputs
        UB(i) = Plant.Generator(i).Size;
        if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
            stor(end+1) = i;
        elseif strcmp(Plant.Generator(i).Source,'Renewable')
            renews(end+1) = i;
        elseif ~strcmp(Plant.Generator(i).Type,'Utility')%if its not storage, renewable, or utility, it is dispatchable
            gens(end+1) = i;
        end
    end
    ngens = length(gens);%number of dispatchable generators
    %sort demand into vectors
    for j = 1:1:length(Out)
        GenDemand.(Out{j}) = zeros(ndisps*24,1);
    end
    %load onoff and generator dispatch
    for i = 1:1:length(gens)
        for t = 1:1:ndisps
            OnOff((t-1)*24+1:t*24,gens(i)) = (Plant.Predicted.GenDisp(2:end,gens(i),t)~=0);
            GenDisp((t-1)*25+1:t*25,:) = Plant.Predicted.GenDisp(:,:,t);
        end
    end
    %load demand, scale costs, renewable generation, and fit A dispatch
    %rate plan is the summer winter rates converted into a 1x24 vector of
    %rates at all times of day
    allrates = Plant.Generator(1).VariableStruct.SumRates(:,1);
    RatePlan = allrates(Plant.Generator(1).VariableStruct.SumRateTable(2,:));
    for t = 1:1:ndisps
        for j = 1:1:length(Out)
                GenDemand.(Out{j})((t-1)*24+1:t*24) = Plant.Data.Demand.(Out{j})((t-1)*4+1:4:(t+23)*4);
        end
        scaleCost((t-1)*24+1:t*24,:) = RatePlan*ones(1,nG);%Time*ones(1,nG);
        RenewE((t-1)*25+1:t*25) = sum(Plant.Predicted.GenDisp(:,renews,t),2);
        FirstDisp((t-1)*25+1:t*25,:) = Plant.Predicted.GenDispcQP(:,:,t);
    end
    
    
    %put all inputs into input matrix
    nfeatures = 3*ngens+1+length(stor)+length(Out);
    lockInputs = zeros(nS*ndisps, nfeatures);
    for i = 1:1:length(Out)
        if strcmp(Out,'E')%power demand is the total demand minus any must take renewable generation
            lockInputs(:,i) = (GenDemand.(Out{i})-RenewE)./(max(GenDemand.(Out{i})));
        else
            lockInputs(:,i) = GenDemand.(Out{i})./(max(GenDemand.(Out{i})));
        end
    end
    col = length(Out)+1;%collumn for input
    %add in cost terms
    lockInputs(:,col) = (scaleCost(:,1).*dt)./max(Plant.Generator(1).VariableStruct.SumRates(:,1)).*dt;%cost of utility electrical
    %add in power from storage given fit A
    for i = 1:1:ndisps
        lockInputs((i-1)*24+1:i*24,col+1:col+length(stor)) = (FirstDisp((i-1)*25+1:i*25-1,stor)-FirstDisp((i-1)*25+2:i*25,stor)+ones(nS,1)*UB(stor))./(ones(nS,1)*(max(FirstDisp(1:end-1,stor)-FirstDisp(2:end,stor))+UB(stor)));%SOC of battery at previous step
    end
    col = col+length(stor)+1;
    for i = 1:1:ngens
        j = gens(i);
        j1 = Plant.OneStep.organize{j};
        if ~isempty(j1)
            for t = 1:1:ndisps
                %lockInputs((t-1)*24+1:t*24,col) = min(FirstDisp((t-1)*25+1:t*25-1,j)+Plant.Generator(j).QPform.Ramp.b(1),UB(j))./UB(j);%upper bound
                %lockInputs((t-1)*24+1:t*24,col+1) = max(FirstDisp((t-1)*25+1:t*25-1,j)-Plant.Generator(j).QPform.Ramp.b(1),LB(j))./UB(j);%lower bound
                lockInputs((t-1)*24+1:t*24,col) = FirstDisp((t-1)*25+2:t*25,j)./UB(j);%estimated dispatch from FitA at this timestep
                lockInputs((t-1)*24+1:t*24,col+1) = (FirstDisp((t-1)*25+1:t*25-1,j)-LB(j))./UB(j);%estimated dispatch from FitA at the previous timestep
                lockInputs((t-1)*24+1:t*24,col+2) = [FirstDisp((t-1)*25+3:t*25,j);FirstDisp((t-1)*25+1,j)]./UB(j);%estimated dispatch from FitA at the next timestep
%                 lockInputs((t-1)*24+1:t*24,col+3) = [FirstDisp((t-1)*25+4:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+2,j)]./UB(j);%estimated dispatch from FitA at t2
%                 lockInputs((t-1)*24+1:t*24,col+4) = [FirstDisp((t-1)*25+5:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+3,j)]./UB(j);%estimated dispatch from FitA at t3
%                 lockInputs((t-1)*24+1:t*24,col+5) = [FirstDisp((t-1)*25+6:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+4,j)]./UB(j);%estimated dispatch from FitA at t4
%                 lockInputs((t-1)*24+1:t*24,col+6) = [FirstDisp((t-1)*25+7:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+5,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+7) = [FirstDisp((t-1)*25+8:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+6,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+8) = [FirstDisp((t-1)*25+9:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+7,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+9) = [FirstDisp((t-1)*25+10:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+8,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+10) = [FirstDisp((t-1)*25+11:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+9,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+11) = [FirstDisp((t-1)*25+12:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+10,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+12) = [FirstDisp((t-1)*25+13:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+11,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+13) = [FirstDisp((t-1)*25+14:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+12,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+14) = [FirstDisp((t-1)*25+15:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+13,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+15) = [FirstDisp((t-1)*25+16:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+14,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+16) = [FirstDisp((t-1)*25+17:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+15,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+17) = [FirstDisp((t-1)*25+18:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+16,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+18) = [FirstDisp((t-1)*25+19:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+17,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+19) = [FirstDisp((t-1)*25+20:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+18,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+20) = [FirstDisp((t-1)*25+21:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+19,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+21) = [FirstDisp((t-1)*25+22:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+20,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+22) = [FirstDisp((t-1)*25+23:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+21,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+23) = [FirstDisp((t-1)*25+24:t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+22,j)]./UB(j);%estimated dispatch from FitA at t5
%                 lockInputs((t-1)*24+1:t*24,col+24) = [FirstDisp(t*25,j);FirstDisp((t-1)*25+1:(t-1)*25+23,j)]./UB(j);%estimated dispatch from FitA at t5               
            end
            col = col+3;
        end
    end

    %prevent Nans
    lockInputs(isnan(lockInputs)) = 0;
    lockInputs = lockInputs(1:8712*nS,:);
    OnOff = OnOff(1:8712*nS,:);
    %convert from 0 to 1 scale to -1 to 1 scale
    lockInputs = lockInputs*2-1;
    %% sort instances into training and testing
    %you want 70% of all instances for training, and 30% for testing, but
    %you want an even spread from all hours and all days, so use 7 as
    %dividing factor
    
    %take out every seventh entry for testing
    testi = [1:9:length(OnOff(:,1))];
    %create vector of everything that isn't pulled out for testing
    traini = zeros(length(OnOff(:,1))-length(testi),1);
    j =1;
    for i = 1:1:length(OnOff(:,1))
        if rem(i,9) ~= 0
            traini(j) = i;
            j = j+1;
        end
    end
    lockInputs_train = lockInputs(traini,:);
    OnOff_train = OnOff(traini,:);
    lockInputs_test = lockInputs(testi,:);
    OnOff_test = OnOff(testi,:);
    
    %random order so that order doesn't effect output
    newOrder_train = randperm(length(traini));
    lockInputs_train = lockInputs_train(newOrder_train,:).*100;
    OnOff_train = OnOff_train(newOrder_train,:);
    
    newOrder_test = randperm(length(testi));
    lockInputs_test = lockInputs_test(newOrder_test,:).*100;
    OnOff_test = OnOff_test(newOrder_test,:);
    

    %% test trees first to find out which inputs are valuable
    Tree = [];
    TrimmedT = [];
    %save inputs and solutions into csv file for weka testing
    %each generator needs a different tree
%     for i = 1:1:length(gens)
%         sols_train = (OnOff_train(:,gens(i))==1);
%         sols_test = (OnOff_test(:,gens(i))==1);
%         trainfilename = strcat('locktrain_',num2str(gens(i)),'.csv');
%         testfilename = strcat('locktest_',num2str(gens(i)),'.csv');
%         csvwrite(trainfilename,[lockInputs_train,sols_train]);
%         csvwrite(testfilename,[lockInputs_test,sols_test]);
%     end

    %% use libsvm to find polynomial kernel with degree of highest accuracy
    percOn = [];
    %test accuracy for each component at each degree
%     degs = [0,1,2,3,4,5];
%     %preallocate accuacy measurements
%     AccTrain = zeros(length(degs), length(gens));
%     AccValid = zeros(length(degs), length(gens));
%     nSVMs = zeros(length(degs), length(gens));%number of support vector machines
%     percOn = zeros(length(gens),1);
%     %test for each component
%     for i = 1:1:length(gens)
%         sols = OnOff_train(:,gens(i));
%         valid = OnOff_test(:,gens(i));
%         %convert to doubles
%         solsd = double(sols);
%         validd = double(valid);
%         %record percent of examples with differing outputs
%         percOn(i) = (nnz(sols)/length(sols))*100;
%         %convert to 1, -1 instead of 0 1
%         solsd(solsd==0) = -1;
%         validd(validd==0) = -1;
%         %if the training set only contains one class
%         if nnz(sols) == length(sols) || nnz(sols) == 0
%             AccTrain(:,i) = 1; 
%             AccValid(:,i) = nan;
%         else
%             %test at each degree
%             for j = 1:1:length(degs)
%                 %setup the model
%                 if degs(j) == 0%if you are training a linear kernel function
%                     trainparams = '-t 0 -s 2 -q';
%                 else%if you are training a polynomial kernel function
%                     trainparams = sprintf('-t 1 -d %d -s 2 -q',degs(j));
%                 end
%                 
%                 %create the model
%                 model = svmtrain(lockInputs_test, validd, trainparams);
%                 %record number of support vectors
%                 nSVMs(j,i) = length(model.sv_indices);%m.get_nr_sv();
%                 %record training accuracy
%                 [~,acci,~]= svmpredict(validd,lockInputs_test, model,'-q');%[~, AccTrain(j,i), ~] = svm_predict(subtrainClass, subtraindata, m);
%                 AccTrain(j,i) = acci(1);
%                 %record test accuracy
%                 [~, acci, ~] = svmpredict(solsd,lockInputs_train, model, '-q');%svm_predict(validateClass, validatedata, m);
%                 AccValid(j,i) = acci(1);
%             end
%         end
%     end
    
    %% check neural network configuration
    
    %chec for all components at all depths in question
    depth = [1:5];%this will become the number of hidden layers: if depth = 1 the network will have two layers according to the toolbox
    %initialize accuracy for training and validation
    AccTrain = zeros(length(depth), length(gens));
    AccValid = zeros(length(depth), length(gens));
    for i = 1:1:length(depth)
        for j = 1:1:length(gens)%one network per generator
            sols_train = OnOff_train(:,gens(j))';
            sols_test = OnOff_test(:,gens(j))';
            %use matlab neural network library
            [ANN, accTrain, accTest] = toolboxANNtest(depth(i), lockInputs_test', lockInputs_train', sols_test, sols_train);
            AccTrain(i,j) = accTrain;
            AccValid(i,j) = accTest;
        end
    end
    
end

