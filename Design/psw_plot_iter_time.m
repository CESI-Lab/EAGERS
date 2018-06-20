function stop = psw_plot_iter_time(optimValues, state)
%PSW_PLOT_ITER_TIME

% persistent variables
persistent startTime

% parameters
pointSize = 15;
color = 'blue';

% optimization shouldn't stop
stop = false;

switch state
    case 'init'
        startTime = doInit();
        startTime = doIter(optimValues, startTime, pointSize, color);
    case 'iter'
        startTime = doIter(optimValues, startTime, pointSize, color);
    case 'done'
        hold off
end

end


%% Local functions

function startTime = doInit()
% start timer
startTime = tic;
% plot setup
title('Iteration Time')
xlabel('Iteration Number')
ylabel('Time for Iteration Completion [min]')
hold on
end

function newStartTime = doIter(optimValues, startTime, pointSize, ...
    color)
x = optimValues.iteration;
y = toc(startTime) / 60;
newStartTime = tic;
scatter(x, y, pointSize, color)
end
