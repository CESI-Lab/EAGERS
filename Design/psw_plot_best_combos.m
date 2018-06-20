function stop = psw_plot_best_combos(optimValues, state)
%PSW_PLOT_BEST_COMBOS

% persistent variables
persistent xComb yComb lastBest

% handle inputs
nDim = length(optimValues.bestx);

% parameters
pointSize = 15;
colorScheme = hsv(nDim);

% optimization shouldn't stop
stop = false;

switch state
    case 'init'
        [xComb, yComb, lastBest] = doInit(colorScheme, nDim);
        [xComb, yComb, lastBest] = doIter(...
            xComb, yComb, lastBest, optimValues, pointSize);
    case 'iter'
        [xComb, yComb, lastBest] = doIter(...
            xComb, yComb, lastBest, optimValues, pointSize);
    case 'done'
        hold off
end

end


%% Local functions

function [x, y, lastComb] = doInit(colorOrder, nDim)
% calculate outputs
x = inf(200*nDim, 1);
y = inf(200*nDim, nDim);
lastComb = NaN;
% plot setup
set(gca, 'ColorOrder', colorOrder)
title('Best Combinations')
xlabel('Iteration Number')
ylabel('State Values')
hold on
end

function [x, y, lastBest] = doIter(x, y, lastBest, ...
    optimValues, pointSize)
% quit if there has been no improvement
thisBest = optimValues.bestx;
if thisBest == lastBest
    return
else
    lastBest = thisBest;
end
% get next index to assign to
iNonInf = find(~isinf(x));
if isempty(iNonInf)
    iNext = 1;
else
    iNext = iNonInf(end) + 1;
end
% assign values
iter = optimValues.iteration;
x(iNext) = iter;
y(iNext,:) = optimValues.bestx;
% plot
h = plot(x, y, 'o');
set(h, 'MarkerSize', pointSize/4)
end
