function [eVented,xL,lb,ub] = vented_energy(Gen,subnet,vent,xL,lb,ub,param)
%add states at each node to allow puroseful energy discharge (venting heat or cooling)
if strcmp(param,'H') 
    net = 'DistrictHeat';
elseif strcmp(param,'C')
    net = 'DistrictCool';
end
if any(strcmp(net,fieldnames(subnet))) && vent == 1
    %%find maximum heat production possible
    maxGen = 0;
    nG = length(Gen);
    for i = 1:1:nG
        if isfield(Gen(i).QPform.output,param)
            if ~isempty(strfind(Gen(i).Type,'Storage'))
                maxGen = maxGen + Gen(i).QPform.Stor.PeakDisch;
            else
                if strcmp(param,'H') && isfield(Gen(i).QPform.output,'E')
                    maxGen = maxGen + max(Gen(i).Size*Gen(i).Output.Capacity./Gen(i).Output.Electricity.*Gen(i).Output.Heat);
                elseif Gen(i).QPform.output.(param)(end)>0
                    maxGen = maxGen + max(max(Gen(i).QPform.output.(param)))*Gen(i).Size;
                end
            end
        end
    end
    %%assume heat can be lost any any node in the network that has a device producing heat
    n = length(subnet.(net).nodes);
    ventNodes = 0;
    eVented =zeros(n,1); %matrix for the state associated with venting heat at each district heating node, at each time step
    for i = 1:1:n
        genI = subnet.(net).Equipment{i};%%identify generators at this node
        for j = 1:1:length(genI)
            if isfield(Gen(genI(j)).QPform.output,'H')
                ventNodes = ventNodes + 1; % Add single state for heat that is ventd to make energy equality true
                eVented(1,i) = (xL+ventNodes);
                break
            end
        end
    end
    xL = xL+ventNodes;
    lb(end+1:end+ventNodes,1) = 0;
    ub(end+1:end+ventNodes,1) = maxGen;
elseif any(strcmp(net,fieldnames(subnet)))
    eVented = zeros(length(subnet.(net).nodes),1);
else
    eVented = [];
end
end%Ends function vented_energy