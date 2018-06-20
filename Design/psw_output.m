function stop = psw_output(optimValues, state)
%PSW_OUTPUT Output function for particleswarm optimization.

% start time persistence for time tracking
persistent startTime

% optimization shouldn't stop
stop = false;

% time tracking
switch state
    case 'init'
        startTime = tic;
    case 'iter'
        printTime(optimValues, toc(startTime))
    case 'done'
        printTime(optimValues, toc(startTime))
end

end


%% Local functions

function printTime(optimValues, totalSeconds)
%PRINTTIME

seconds = mod(totalSeconds, 60);
minutes = mod(floor(totalSeconds / 60), 60);
hours = floor(totalSeconds / 3600);
fprintf('  Iteration %03i  --  %2.0f : %2.0f : %4.1f    elapsed.\n', ...
    optimValues.iteration, hours, minutes, seconds)

end
