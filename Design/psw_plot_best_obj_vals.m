function stop = psw_plot_best_obj_vals(optimValues, state)
%PSW_PLOT_BEST_OBJ_VALS

% parameters
pointSize = 15;
color = 'blue';

% optimization shouldn't stop
stop = false;

switch state
    case 'init'
        doInit()
        doIter(optimValues, pointSize, color)
    case 'iter'
        doIter(optimValues, pointSize, color)
    case 'done'
        hold off
end

end


%% Local functions

function doInit()
% plot setup
title('Best Objective Function Value')
xlabel('Iteration Number')
ylabel('Objective Function Value')
hold on
end

function doIter(optimValues, pointSize, color)
% determine whether objective function has returned any Inf's
if any(any(isinf(optimValues.swarmfvals)))
    color = 'red';
end
% assign values
x = optimValues.iteration;
y = optimValues.bestfval;
% plot
scatter(x, y, pointSize, color)
end
