function time = build_time_vector(options)
%% set up dt vector of time interval length
if strcmp(options.tspacing,'constant') == 1
    n_s = options.Horizon/options.Resolution; % # of time intervals
    time = linspace(options.Resolution,options.Horizon,n_s)';
%% if user manually specified timesteps
elseif strcmp(options.tspacing, 'manual')
    time = options.manualT.*options.Horizon; %the user was told to specify steps as a portion of the horizon
    if options.manualT(end)~=1 %if they did not specify up to one horizon make the last entry one horizon
        time(end+1) = options.Horizon;
    end
    round_time = round(time/(options.Tmpc/3600));%make sure that each step can contain a whole number of MPC steps
    time = round_time*(options.Tmpc/3600);
    [m,n] = size(time);
    if m==1 && n>1
        time = time';
    end
%% if user specified logarithmic timesteps, then have more steps in the beggining, placing less certainty on the distant predictions
elseif strcmp(options.tspacing, 'logarithm')
    %the first step must be longer than the time for tmpc
    %made the first step a 5 min interval
    n1 = round(options.Horizon);
    if options.Topt/3600>1/12 %if the online loop is greater than 5 minutes, then make the shortest interval one onlineloop
        const = 1/(options.Topt/3600);
    else
        const = 12; %if the online loop runs more frequently than 5 min, then make the shortest interval 5 min
    end
    t = exp(1:n1)'/(const*exp(1));%divide by 12, because you want the first step to be 5 min
    time = t(t<=24);
    round_time = round(time/(options.Tmpc/3600));%make sure that each step can contain a whole number of MPC steps
    time = round_time*(options.Tmpc/3600);
end
end%Ends function build_time_vector