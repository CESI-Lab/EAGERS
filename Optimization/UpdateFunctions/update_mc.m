function marginal = update_mc(gen,dispatch,scale_cost,dt,v_h)
%%calculate the marginal cost for any networks it can be computed for
n_g = length(gen);
marginal = [];
ac_dc =[];
s = {};
for i = 1:1:n_g
    switch gen(i).Type
        case 'Utility'
            if ~isempty(gen(i).QPform.output)
                f = fieldnames(gen(i).QPform.output);
                if ~any(strcmp(f{1},s)) 
                    s(end+1) = f(1);
                end
            end
        case {'Electric Storage';'Electric Generator';}
             if isfield(gen(i).Output,'Electricity') && ~any(strcmp('E',s))   
                 s(end+1) = {'E'};
             elseif isfield(gen(i).Output,'DirectCurrent') && ~any(strcmp('DC',s))   
                 s(end+1) = {'DC'};
             end
        case 'CHP Generator'
            if isfield(gen(i).Output,'Electricity') && ~any(strcmp('E',s))   
                s(end+1) = {'E'};
            elseif isfield(gen(i).Output,'DirectCurrent') && ~any(strcmp('DC',s))   
                s(end+1) = {'DC'};
            end
            if ~any(strcmp('H',s))   
                s(end+1) = {'H'};
            end
        case 'Heater'
            if ~any(strcmp('H',s))   
                s(end+1) = {'H'};
            end
        case 'Chiller'
            if ~any(strcmp('C',s)) 
                s(end+1) = {'C'};
            end
        case 'Thermal Storage'
            if isfield(gen(i).QPform.output,'H') && ~any(strcmp('H',s)) 
                s(end+1) = {'H'};
            elseif isfield(gen(i).QPform.output,'C') && ~any(strcmp('C',s)) 
                s(end+1) = {'C'};
            end
        case 'Hydrogen Storage'
            if ~any(strcmp('Hy',s))   
                s(end+1) = {'Hy'};
            end
        case 'Hydro Storage'
            if ~any(strcmp('W',s))   
                s(end+1) = {'W'};
            end
        case 'Cooling Tower'
            if ~any(strcmp('CW',s))   
                s(end+1) = {'CW'};
            end
        case 'AC_DC'
            ac_dc = gen(i);
            if ~any(strcmp('E',s))   
                s(end+1) = {'E'};
            end
            if ~any(strcmp('DC',s))   
                s(end+1) = {'DC'};
            end
    end
end

%determine if it needs to split the marginal cost of CHP generators between heating and electric
if any(strcmp('E',s)) || any(strcmp('DC',s)) || any(strcmp('H',s))
    chp = [];
    for i = 1:1:n_g
        if strcmp(gen(i).Type,'CHP Generator')
            chp(end+1) = i;
        end
    end
    if ~isempty(chp) && isempty(v_h)
        v_h = true;
    end
end

%%Electricity or Direct current marginal cost
if any(strcmp('E',s)) || any(strcmp('DC',s))
    include = [];
    dc_to_ac = ones(1,n_g);
    ac_to_dc = ones(1,n_g);
    for i = 1:1:n_g
        if gen(i).Enabled && (strcmp(gen(i).Type,'CHP Generator') || strcmp(gen(i).Type,'Electric Generator') || strcmp(gen(i).Type,'Hydrogen Generator'))
            include(end+1) = i;
            if isfield(gen(i).Output,'DirectCurrent') && ~isempty(ac_dc)
                dc_to_ac(i) = 1/ac_dc.VariableStruct.DC_to_AC_eff;
            end
            if isfield(gen(i).Output,'Electricity') && ~isempty(ac_dc)
                ac_to_dc(i) = 1/ac_dc.VariableStruct.AC_to_DC_eff;
            end
        end
    end
    scale = scale_cost;
    if v_h
        scale(:,chp) = scale(:,chp)*.7;
    end
    if isempty(dispatch) || length(dispatch(:,1))==1
        if any(strcmp('E',s))
            marginal.E = m_loop(gen,include,scale.*dc_to_ac,dispatch,'E',[]);
        end
        if any(strcmp('DC',s))
            marginal.DC = m_loop(gen,include,scale.*ac_to_dc,dispatch,'DC',[]);
        end
    else
        if any(strcmp('E',s))
            for t = 1:1:length(dt)
                scale(t,:) = scale(t,:).*dc_to_ac;
            end
            marginal.E =  min_max_cost(gen,include,scale,dispatch,'E',[],dt);
        end
        if any(strcmp('DC',s))
            for t = 1:1:length(dt)
                scale(t,:) = scale(t,:).*ac_to_dc;
            end
            marginal.DC =  min_max_cost(gen,include,scale,dispatch,'E',[],dt);
        end
    end
end


if any(strcmp('H',s))
    include = [];
    for i = 1:1:n_g
        if gen(i).Enabled && (strcmp(gen(i).Type,'CHP Generator') || strcmp(gen(i).Type,'Heater'))
            include(end+1) = i;
        end
    end
    scale = scale_cost;
    if v_h
        scale(:,chp) = scale(:,chp)*.3;
    else
        scale(:,chp) = scale(:,chp)*.1;
    end
    if isempty(dispatch) || length(dispatch(:,1))==1
        marginal.H = m_loop(gen,include,scale,dispatch,'H',[]);
    else
        marginal.H = min_max_cost(gen,include,scale,dispatch,'H',[],dt);
    end
end

if any(strcmp('C',s))%%chillers are unique because the cost in QP.f is zero
    include = [];
    for i = 1:1:n_g
        if gen(i).Enabled && strcmp(gen(i).Type,'Chiller')
            include(end+1) = i;
        end
    end
    if isempty(dispatch) || length(dispatch(:,1))==1
        marginal.C = m_loop(gen,include,scale_cost,dispatch,'C',marginal);
    else
        marginal.C = min_max_cost(gen,include,scale_cost,dispatch,'C',marginal,dt);
    end
end

if any(strcmp('W',s))
    if isempty(dispatch) || length(dispatch(:,1))==1
        marginal.W = 1;
    else
        marginal.W.Min = 0;
        marginal.W.Max = 1;
    end
end

if any(strcmp('Hy',s)) 
    marginal.Hy = marginal.E;
%     if isempty(Dispatch) || length(Dispatch(:,1))==1
%         marginal.Hy = 1;
%     else
%         marginal.Hy.Min = 0;
%         marginal.Hy.Max = 1;
%     end
end
end%Ends function update_mc

function mc = m_loop(gen,include,scale,dispatch,s,marginal)
n_g = length(gen);
margin_cost = [];
for j = 1:1:length(include)
    i = include(j);
    if isempty(dispatch)
        margin_cost(end+1) = m_cost(gen(i).QPform,[],scale(i),marginal);
    else
        margin_cost(end+1) = m_cost(gen(i).QPform,dispatch(1,i),scale(i),marginal);
    end
end
mc = nan;
if ~isempty(dispatch)
    for j = 1:1:length(include)
        i = include(j);
        states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
        if dispatch(1,i)>gen(i).QPform.(states{1}).lb(end)%on above lower bound
            mc = max(margin_cost(j),mc);
        end
    end
end
if isempty(dispatch) || isnan(mc)%  || mc == 0 
    mc = min(margin_cost);
end
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Utility') && isfield(gen(i).QPform.output,s)
        if isempty(mc) || gen(i).QPform.X.f*scale(i)<mc
            mc = gen(i).QPform.X.f*scale(i);
        end
    end
end
end%Ends function m_loop

function m_c = m_cost(gen,set,scale,marginal)
states = gen.states(1:nnz(~cellfun('isempty',gen.states(:,end))),end);
j = 1;
m_c = gen.(states{1}).f(end)*scale;
if ~isempty(set)     
    while j<length(states) && gen.(states{j}).ub(end)<=set
        set = set - gen.(states{j}).ub(end);
        j = j+1;
    end
    m_c = (gen.(states{j}).f(end) + set*gen.(states{j}).H(end))*scale;
end
if m_c == 0%things like chillers where input cost is output of a different generator
    f = fieldnames(gen.output);
    for s = 1:1:length(f)
        if length(gen.output.(f{s})(1,:))>1 && all(all(gen.output.(f{s})<=0))
            m_c = marginal.(f{s})*(-gen.output.(f{s})(j,end));
        end
    end
end
if m_c<0
    j = 1;
    while j<length(states) && m_c<0
        j = j+1;
        m_c = gen.(states{j}).f(end)*scale;
    end
    if m_c<0
        m_c = 1e-3;
    end
end
end%Ends function m_cost

function mc = min_max_cost(gen,include,scale,dispatch,s,marginal,dt)
n_g = length(gen);
n_s = length(dt);
min_c = zeros(n_s,length(include));
max_c = zeros(n_s,length(include));
charge_index = false(n_s,1);
for k = 1:1:length(include)
    i = include(k);
    states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
    min_c(:,k) = gen(i).QPform.(states{1}).f(end)*scale(:,i);
    max_c(:,k) = (gen(i).QPform.(states{end}).f(end) + gen(i).QPform.(states{end}).ub(end)*gen(i).QPform.(states{end}).H(end))*scale(:,i);
    if all(max_c(:,k)== 0) %handles chillers or other things that don't have a direct cost in the cost function
        f = fieldnames(gen(i).QPform.output);
        for s = 1:1:length(f)
            if all(gen(i).QPform.output.(f{s})(:,end)<=0)
                min_c(:,k) = marginal.(f{s}).Min*min(-gen(i).QPform.output.(f{s})(:,end));
                max_c(:,k) = marginal.(f{s}).Max*max(-gen(i).QPform.output.(f{s})(:,end));
            end
        end
    end
    if any(min_c(:,k)<0)
        j = 1;
        while j<length(states) && any(min_c(:,k)<0)
            j = j+1;
            min_c(:,k) = gen(i).QPform.(states{j}).f(end)*scale(:,i);
        end
        if any(min_c(:,k)<0)
            min_c(:,k) = 1e-3;
        end
    end
%     max_c(:,k) = min(3*min_c(:,k),max_c(:,k));
end
for i = 1:1:n_g
    if isfield(gen(i).QPform,'Stor') && isfield(gen(i).QPform.output,s)
        charge_index = max(charge_index,(dispatch(2:end,i)-dispatch(1:end-1,i))>0);
    end
end
if any(any((dispatch(2:end,include)>0)))
    min_c(dispatch(2:end,include)==0) = nan;
    max_c(dispatch(2:end,include)==0) = nan;
end
gen_on_during_charge = false(n_s,length(include));
charge_index2 = nonzeros((1:n_s)'.*charge_index)+1;
gen_on_during_charge(charge_index,:) = dispatch(charge_index2,include)>0;
if any(charge_index>0) && any(any(gen_on_during_charge))
    max_c(gen_on_during_charge==1) = 1.25*max_c(gen_on_during_charge==1); %if the storage is charging, then its margin cost can be greater than the cost of the generators that are on
end
min_on_t = min(min_c,[],2); 
max_on_t = max(max_c,[],2);
timedivide_min = max(1,sum(dt(min_on_t~=0)));
timedivide_max = sum(dt(max_on_t>0));
min_on_t(isnan(min_on_t))=0;
max_on_t(isnan(max_on_t))=0;
mc.Min = sum(min_on_t(min_on_t~=0).*dt(min_on_t~=0))/timedivide_min;%make the cost proportional to amount of time at that cost
mc.Max = sum(max_on_t(max_on_t>0).*dt(max_on_t>0))/timedivide_max;

for i = 1:1:n_g
    if strcmp(gen(i).Type,'Utility') && isfield(gen(i).QPform.output,s)
        if isfield(gen(i).QPform,'Y') %has sellback
            mc.Min = min(mc.Min,-gen(i).QPform.Y.f*min(scale(:,i)));
        else
            mc.Min = min(mc.Min,gen(i).QPform.X.f*min(scale(:,i)));
        end
        mc.Max = min(mc.Max,gen(i).QPform.X.f*max(scale(:,i)));
    end
end
end%Ends function min_max_cost