function [dx_dt,flag] = ss_response(gen,handles)
% Flag values:
%   0 -- Standard operation.
%   1 -- error calculating response curve.
switch gen.Type
        case {'CHP Generator';'Electric Generator';'Hydrogen Generator';} 
            names = {'Electricity';};
            if isfield(gen.VariableStruct.Startup,'Electricity')
                lower_bound = gen.VariableStruct.Startup.Electricity(end);
            elseif isfield(gen.VariableStruct.Startup,'DirectCurrent')
                lower_bound = gen.VariableStruct.Startup.DirectCurrent(end);
            end
        %     if strcmp(gen.Type,'CHP Generator')
        %        Names = {'Electricity';'Heat';}; 
        %     end
    case'Electrolyzer'
        names = {'Hydrogen';};
        lower_bound = gen.VariableStruct.Startup.Hydrogen(end);
    case 'Heater'
        names = {'Heat'};
        lower_bound = gen.VariableStruct.Startup.Heat(end);
    case 'Chiller'
        names = {'Cooling'};
        lower_bound = gen.VariableStruct.Startup.Cooling(end);
    case 'Cooling Tower'
        lower_bound = gen.VariableStruct.Startup.heat_reject(end);
end
if isfield(gen.VariableStruct,'dX_dt')
    t_peak = (gen.Size-lower_bound)/gen.VariableStruct.dX_dt*3600;
else 
    t_peak = 3600;
end
if t_peak>4*60
    dt = 60;
elseif t_peak>3.6e4
    dt = 3600;
else
    dt = 1;
end
%convert continuous state-space to discrete time state-space
a = expm(gen.VariableStruct.StateSpace.A.*dt);
b = (gen.VariableStruct.StateSpace.A)\(expm(gen.VariableStruct.StateSpace.A.*dt)-eye(2))*gen.VariableStruct.StateSpace.B;
c = gen.VariableStruct.StateSpace.C;
d = gen.VariableStruct.StateSpace.D;
time = 5*t_peak;
dx_dt = [];
count = 0;
flag = 0;
while isempty(dx_dt)
    nS = round(time/dt)+1;
    u = gen.Size*linspace(1,1,nS);
    [n,n2] = size(gen.VariableStruct.StateSpace.C);
    x0 = zeros(n2,1);
    for i = 1:1:n
        x0(find(gen.VariableStruct.StateSpace.C(i,:),1,'first'))=lower_bound(i);
    end
    [y,t] = ss_sim(a,b,c,d,x0,u,dt);
    if any((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1)))) && ~any(abs((y(end-9:end,1)-y(end,1))/y(end,1))>1e-3)
        nR = min(nonzeros((1:1:nS+1)'.*((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1))))));
        t_rise = interp1(y(nR-1:nR,1)-y(1,1),t(nR-1:nR),.95*(u(1)-y(1,1)))/3600; %rise time in hours
        dx_dt = ((gen.Size-lower_bound).*(0.95)./t_rise);
    else
        time = time*5;
    end
    count = count+1;
    if count>10
        disp('error in ss_response');
        flag = 1;
        dx_dt = 1;
%         nR = min(nonzeros((1:1:nS+1)'.*((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1))))));
%         t_rise = interp1(y(nR-1:nR,1)-y(1,1),t(nR-1:nR),.95*(u(1)-y(1,1)))/3600; %rise time in hours
%         dx_dt = ((gen.Size-lower_bound).*(0.95)./t_rise);
    end
end
if ~isempty(handles)
    h1 = handles.ResponseRate;
    cla(h1)
    if dt==1
        plot(h1,t,y);
        xlabel(h1,'Time (s)')
    elseif dt==60
        plot(h1,t/60,y);
        xlabel(h1,'Time (min)')
    elseif dt==3600
        plot(h1,t/3600,y);
        xlabel(h1,'Time (hr)')
    end
    ylabel(h1,'Power (kW)')
    % set(get(h1,'Ylabel'),'String','Power (kW)')
    legend(h1, names)
end
end%Ends function second_order_response

function [y,time] = ss_sim(a,b,c,d,x0,u,d_t)
n = length(u);
time = (0:d_t:n*d_t)';
x = zeros(length(x0),n+1);
x(:,1) = x0;
y = zeros(n+1,1);
y(1) = c*x0;
for t = 1:1:n
    x(:,t+1) = a*x(:,t) + b*u(t);
    y(t+1) = c*x(:,t+1) + d*u(t);
end
end%Ends function ss_sim
