function dY = RunBlocks(t,Y)
global modelParam IterCount Inlet Outlet TagInf TagFinal Tags SimSettings WaitBar Jcount
Y = Y.*modelParam.Scale;
CompNames = fieldnames(modelParam.Components);
controls = fieldnames(modelParam.Controls);
list = [CompNames;controls;];

if isempty(TagInf) %create TagInf Fields
    TagInf.Time(1) =0;
    TagInf.Pstates = Y(modelParam.Pstates)';
    Jcount = 0;
    for k = 1:1:length(list)
        block = list{k};
        if any(strcmp(controls,block))
            Co = 'Controls';
        else Co = 'Components';
        end
        if isfield(modelParam.(Co).(block),'TagInf')
            f = modelParam.(Co).(block).TagInf;
            for i = 1:1:length(f)
                TagInf.(block).(f{i}) = [];
            end
        end
    end
end
if t==0 || t == TagInf.Time(IterCount)
    Jcount = Jcount+1;
elseif t>TagInf.Time(IterCount)%+1e-5
    IterCount =IterCount+1;
    Jcount = length(Y);
elseif IterCount>1 && t<TagInf.Time(IterCount)
    Jcount = 1;
    while t<TagInf.Time(IterCount-1)
        IterCount = IterCount-1;
    end
    TagInf.Time = TagInf.Time(1:IterCount);
end
TagInf.Time(IterCount,1) = t;

OldInlet = [];%Inlet; %initial condition inlets
%% Update pressure states & inlets
for i = 1:1:length(modelParam.Pstates)
    Pnew = Y(modelParam.Pstates(i));
    if ~isempty(modelParam.Poutlets{i})
        block = modelParam.Poutlets{i}{1};
        port = modelParam.Poutlets{i}{2};
        Outlet.(block).(port) = Pnew; %update outlets based on calculated pressure
    end
end

%%avoid exceeding saturation
Yoriginal = Y;
satLow = Y< modelParam.LowerBound;
Y(satLow) = modelParam.LowerBound(satLow);
satHigh = Y>modelParam.UpperBound;
Y(satHigh) = modelParam.UpperBound(satHigh);
%% run blocks to find outlet conditions & connect inlets
nComp = length(list);
blockSteady = false(nComp,1);
while any(~blockSteady)
    for blockCount = 1:1:nComp 
        block = list{blockCount};
        Inlet.(block) = RefreshInlet(block,t);
        if ~isfield(OldInlet,block)
            OldInlet.(block) = [];
        end
        if HasInletChanged(Inlet.(block),OldInlet.(block))
            OldInlet.(block) = Inlet.(block); %next time only run if the inlets have changed from now.
            if any(strcmp(controls,block))
                Co = 'Controls';
            else Co = 'Components';
            end
            Outlet.(block) = feval(modelParam.(Co).(block).type,t,Y(modelParam.(Co).(block).States),Inlet.(block),modelParam.(Co).(block),'Outlet');
            blockSteady(blockCount) = false;
        else
            blockSteady(blockCount) = true;
        end
    end
end

%% run blocks to find dY
dY = 0*Y;
list = [CompNames;controls;];
for k = 1:1:length(list) %run components with states %% record All tags
    block = list{k};
    if any(strcmp(controls,block))
        Co = 'Controls';
    else Co = 'Components';
    end
    if ~isempty(modelParam.(Co).(block).IC)%blocks with states
        dY(modelParam.(Co).(block).States) = feval(modelParam.(Co).(block).type,t,Y(modelParam.(Co).(block).States),Inlet.(block),modelParam.(Co).(block),'dY');
    end
    if isfield(modelParam.(Co).(block),'TagInf')
        tagNames = modelParam.(Co).(block).TagInf;
        for i = 1:1:length(tagNames)
            if isempty(TagInf.(block).(tagNames{i}))
                TagInf.(block).(tagNames{i}) = Tags.(block).(tagNames{i});
            elseif isnumeric(Tags.(block).(tagNames{i}))
                TagInf.(block).(tagNames{i})(IterCount,:)=Tags.(block).(tagNames{i});
            elseif isstruct(Tags.(block).(tagNames{i}))
                f = fieldnames(Tags.(block).(tagNames{i}));
                for j = 1:1:length(f)
                    TagInf.(block).(tagNames{i}).(f{j})(IterCount,:)=Tags.(block).(tagNames{i}).(f{j}); 
                end
            end
        end
    end
    if t== SimSettings.RunTime
        if isfield(modelParam.(Co).(block),'TagFinal')
            tagNames = modelParam.(Co).(block).TagFinal;
            for i = 1:1:length(tagNames)
                if isnumeric(Tags.(block).(tagNames{i}))
                    TagFinal.(block).(tagNames{i})(IterCount,:)=Tags.(block).(tagNames{i});
                else
                    f = fieldnames(Tags.(block).(tagNames{i}));
                    for j = 1:1:length(f)
                        TagFinal.(block).(tagNames{i}).(f{j})(IterCount,:)=Tags.(block).(tagNames{i}).(f{j}); 
                    end
                end
            end
        end
    end
end
TagInf.Pstates(IterCount,:) = Y(modelParam.Pstates)';
dY = dY./modelParam.Scale;
dY = dY+(Y-Yoriginal)./modelParam.Scale; %anti-windup for saturated states

if t>0 && WaitBar.Show == 1 && Jcount==length(Y) && isfield(modelParam,'Scope')
    n = length(modelParam.Scope);
    dt = TagInf.Time(IterCount,1) - TagInf.Time(IterCount-1,1);
    points = nnz(TagInf.Time>(TagInf.Time(IterCount,1)-100*dt));
    for i = 1:1:n
        figure(i)
        tagName = modelParam.Scope{i};
        r = strfind(tagName,'.');
        if isempty(r)
            plot(TagInf.Time(IterCount-points+1:IterCount,1),TagInf.(tagName)(IterCount-points+1:IterCount,:))
        else
            block = tagName(1:r-1);
            tag = tagName(r+1:end);
            plot(TagInf.Time(IterCount-points+1:IterCount,1),TagInf.(block).(tag)(IterCount-points+1:IterCount,:))
        end
        ylabel(tagName);
        xlabel('Time in seconds');
    end
end
%% update waitbar
if t ==0
    Text = 'calculating Jacobian';
    x = Jcount/length(Y);
elseif IterCount>1 && Jcount<length(Y)
    Text = 're-calculating Jacobian';
    x = Jcount/length(Y);
else
    if strcmp(WaitBar.Text,'Running non-linear model with transient')
        x = t/SimSettings.RunTime;
    else x = max(0,(log(t/1e-4))/(log(SimSettings.RunTime/1e-4)));
    end
    Text = {WaitBar.Text;strcat('    Time = ', num2str(t))};
end
if WaitBar.Show == 1
    waitbar(x,WaitBar.Handle,Text);
end
end %ends function RunBlocks

function Change = HasInletChanged(New,Old)
Change = false;
if isempty(Old)
    Change = true;
else
    list2 = fieldnames(New);
    for i = 1:1:length(list2)
        port = list2{i};
        if ~Change  %skip if we already know it has changed
            if isnumeric(New.(port))
                Change = comparePort(Old.(port),New.(port));
            elseif isstruct(New.(port))
                f = fieldnames(New.(port));
                for j = 1:1:length(f)
                    if ~isfield(Old.(port),f{j})
                        Change = true;
                    elseif ~Change %skip if we already know it has changed
                        if isnumeric(New.(port).(f{j}))
                            Change = comparePort(Old.(port).(f{j}),New.(port).(f{j}));
                        else
                            g = fieldnames(New.(port).(f{j}));
                            for k = 1:1:length(g)
                                if ~isfield(Old.(port).(f{j}),g{k})
                                    Change = true;
                                elseif ~Change %skip if we already know it has changed
                                    Change = comparePort(Old.(port).(f{j}).(g{k}),New.(port).(f{j}).(g{k}));
                                end
                            end
                        end
                    end
                end
            elseif iscell(New.(port))
                for j = 1:1:length(New.(port));
                    if New.(port){j}~=Old.(port){j}
                        Change = 1;
                    end
                end
            end
        end
    end
end
end %ends function HasInletChanged

function C = comparePort(Old,New)
tolerance = eps;%1e-14;
if Old ~= 0
    denom = Old;
else denom = 1;
end
Error = (New - Old)/denom;
if  abs(Error) >= tolerance;
    C = true;
else C = false;
end
end %ends function comparePort

function InletBlock = RefreshInlet(block,t)
global modelParam Outlet Tags LinMod
controls = fieldnames(modelParam.Controls);
if any(strcmp(controls,block))
    Co = 'Controls';
else Co = 'Components';
end
list = modelParam.(Co).(block).InletPorts;
for i = 1:1:length(list)
    port = list{i};
    
    if ~isempty(modelParam.(Co).(block).(port).connected)%inlet connected to another outlet
        BlockPort = char(modelParam.(Co).(block).(port).connected);
        if isnumeric(BlockPort);
            InletBlock.(port) = BlockPort;
        else
            r = strfind(BlockPort,'.');
            if ~isempty(r)
                connectedBlock = BlockPort(1:r-1);
                if strcmp(connectedBlock,'Tags')
                    connectedBlock = BlockPort(r(1)+1:r(2)-1);
                    connectedPort = BlockPort(r(2)+1:end);
                    InletBlock.(port) = Tags.(connectedBlock).(connectedPort);
                else
                    connectedPort = BlockPort(r+1:end);
                    InletBlock.(port) = Outlet.(connectedBlock).(connectedPort);
                end
            else
                InletBlock.(port) = feval(BlockPort,t);%finds value of look-up function at time = 0
            end
        end
    else
        InletBlock.(port) = modelParam.(Co).(block).(port).IC;
    end
end
%if linearizing the model, reset the perturbed outlet or tag or lookup function 
if ~isempty(LinMod) && ~isempty(LinMod.Interupt) && strcmp(block,LinMod.Interupt.block)
    if ~isempty(LinMod.Interupt.struct)
        InletBlock.(LinMod.Interupt.port).(LinMod.Interupt.struct)(LinMod.Interupt.index) = LinMod.Interupt.value;
    else
        InletBlock.(LinMod.Interupt.port)(LinMod.Interupt.index) = LinMod.Interupt.value;
    end
end
end %ends function RefreshInlet

function yesNan = IsInletNaN(Inlet)
%% use to help find algebraic loops during initialization
yesNan = false;
list2 = fieldnames(Inlet);
for i = 1:1:length(list2)
    port = list2{i};
    if isnumeric(Inlet.(port))
        yesNan = isnan(Inlet.(port));
    elseif isstruct(Inlet.(port))
        f = fieldnames(Inlet.(port));
        for j = 1:1:length(f)
            if isnumeric(Inlet.(port).(f{j}))
                yesNan =isnan(Inlet.(port).(f{j}));
            else
                g = fieldnames(Inlet.(port).(f{j}));
                for k = 1:1:length(g)
                    yesNan = isnan(Inlet.(port).(f{j}).(g{k}));
                end
            end
        end
    end
end
end %ends function IsInletNaN