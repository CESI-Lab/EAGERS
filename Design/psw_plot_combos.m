function stop = psw_plot_combos(optimValues, state)
%PSW_PLOT_COMBOS

% persistent variables
persistent swarmSize

% handle inputs
nDim = length(optimValues.bestx);

% parameters
pointSize = 15;
colorScheme = hsv(nDim);

% optimization shouldn't stop
stop = false;

switch state
    case 'init'
        swarmSize = doInit(optimValues, colorScheme);
        doIter(swarmSize, optimValues, pointSize);
    case 'iter'
        doIter(swarmSize, optimValues, pointSize);
    case 'done'
        hold off
end

end


%% Local functions

function swarmSize = doInit(optimValues, colorOrder)
% calculate outputs
swarmSize = size(optimValues.swarm, 1);
% plot setup
set(gca, 'ColorOrder', colorOrder)
title('Searched Combinations')
xlabel('Iteration Number')
ylabel('State Values')
hold on
end

function doIter(swarmSize, optimValues, pointSize)
% get next index to assign to
iter = optimValues.iteration;
iLast = iter*swarmSize;
iNext = (iter+1) * swarmSize;
% assign values
x = iter * ones(iNext-iLast,1);
y = optimValues.swarm;
% plot
h = plot(x, y, 'o');
set(h, 'MarkerSize', pointSize/4)
end
