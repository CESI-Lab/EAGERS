function [Net,sqrerror] = trainNetwork(Net,desiredOut, inputs)
%this does forward propagation for a one layer network for a set of
%generators and a demand
%inputs: network, desiredOut: desired network output when using forward
%funnction in the form of a vertical vector, inputs: matrix of inputs of
%size inputlength x number of outputs
%inputs in order: ub, lb, f, H for each generator, demand, $/kWgrid

[sqrerror,dedW,dedb] = finderror(Net,inputs,desiredOut);%find the error and the gradient of the error

%% perturb weights to check that derrordW is found correctly before trying to train using derrordW
% numgrad = zeros(size(Net.Wlayer1));
% perturb = numgrad;
% Winitial = Net.Wlayer1;
% netup = Net;
% netdwn = Net;
% for j = 1:1:length(numgrad(1,:))
%     for i = 1:1:length(numgrad(:,1))
%         perturb(i,j) = .0001*1/Net.nodeconst;
%         netup.Wlayer1 = Winitial+perturb;
%         netdwn.Wlayer1 = Winitial-perturb;
%         [sqrerup,~,~] = finderror(netup, inputs, desiredOut);
%         [sqrerdwn,~,~] = finderror(netdwn, inputs, desiredOut);
%         numgrad(i,j) = (sum(sum(sqrerup-sqrerdwn)))/(2*perturb(i,j));
%         perturb(i,j) = 0;
%     end
% end
% gradientdiff = norm(dedW-numgrad)/norm(dedW+numgrad);
% if gradientdiff>1e-6
%     disp('Warning, gradient function may be incorrect. Training may be inaccurate.')
% end

%% now that you have already determined that derrordW is calculated correctly
%% use BFGS technique to train using derrordW
tolerance = .0001;
% initialize an approximation of the Hessian matrix = d^2f/(dx_i dx_j)
warning('off','all')%prevent print of warning as Hessian gets close to singular
%for i = 1:1:length(dedW(1,:)) %each set for each node output must be trained individually
    iterations = 0;
    laststep = zeros(size(dedW));
    lastbstep = zeros(size(dedb));
    a = 1;
    momentum = .25;%.3 is too high, .1 is too low, .2 does well for test2E_1BS
    while nnz(sqrerror>tolerance)>0 %keep training until you get the desired output
        %find error and relation to weights and biases
        [sqrerror, dedW, dedb] = finderror(Net, inputs, desiredOut);
        iterations = iterations+1;
        step = dedW.*a/100+laststep.*momentum/100;%training step for weights
        bstep = dedb.*a/100+lastbstep.*momentum/100;%training step for bias
        
        Net.Wlayer1 = Net.Wlayer1-step;%minimize error, so go down the slope
        Net.blayer1 = Net.blayer1-bstep;
        
        %check error with new weight an bias
        [sqrerrornew, ~, ~] = finderror(Net, inputs, desiredOut);
        
        %if it gets worse, try the other direction and try a different step size 
        if sum(sum(abs(sqrerrornew)))>=sum(sum(abs(sqrerror))) || nnz(isinf(sqrerrornew))>0%|| nnz(isnan(dedWnew))>0%if the error gets worse or you have reached a flat point
            if abs(a) <1e-12 %try different size steps
                %direction = direction + 1;
                a = 1;
            else
                a = a/10;
            end
            %undo the last change
            Net.Wlayer1 = Net.Wlayer1+step;
            Net.blayer1 = Net.blayer1+bstep;
            laststep = zeros(size(laststep));
            lastbstep = zeros(size(lastbstep));
            
        %if it gets better, keep the change and keep going
        else %if it works 
            laststep = step;
            lastbstep = bstep;
            sqrerror = sqrerrornew;
            %if you are below tolerance, go to the next weight
            if nnz(sqrerrornew>tolerance)==0
                disp('below tolerance');
                break
            end
        end
        
        %if you have hit your maz iterations, stop
        if iterations>1e+4 %|| direction >=2 %|| (iterations>1e+4*length(desiredOut(1,:))/100 && nnz(sqrerrornew1>tolerance)/length(desiredOut(1,:))<0.01)%/100%10^4 iterations per 100 timesteps
            disp('not converging after 10^4 iterations, exiting loop');
            if Net.classify
                sqrerror = sqrerrornew;
                break
            else
                sqrerror = sqrerrornew;
                break
            end
        end
        
    end



function [cost,derrordW,derrordb] = finderror(Net,inputs,desiredOut)
NetOut = forward(Net,inputs);
error = (desiredOut-NetOut);%all errors
cost = error.^2.*0.5;%/length(inputs(:,1)) + (Net.lambda/2)*sum(Net.Wlayer1)'.^2*ones(1,length(desiredOut(1,:)));% normalized(1/2error^2) + 1/2penalty for complex model
%keep model simple using lambda to prevent over fitting
if Net.classify %if it has a sigmoid function
    %use cross error to prevent learning slowdown with sigmoid functions
    derrordW = -2*inputs'*(error.*(NetOut.*(1-NetOut))*Net.nodeconst)/length(desiredOut(:,1));
    derrordb = -2*sum(error.*(NetOut.*(1-NetOut))*Net.nodeconst)/length(desiredOut(:,1));
else%if no activation function
    derrordW = (-error*inputs)';% + Net.lambda*Net.Wlayer1;%no activation function so this is just the error
    derrordb = -1/length(error(1,:))*sum(error,2)';
end

