function optimize_component_size(gen_i,min_size,max_size)
%this function re-sizes a single component to reduce the NPC of the
%building operation
global testSystems SYSINDEX  TestData
handles = guihandles;
n = 5;
gen_i_size = linspace(min_size,max_size,n);

if get(handles.DesignDay,'Value') == 1
    design_day = true;
else
    design_day = false;
end
years = str2double(get(handles.NPC_Years,'String'));%choose years in GUI
% find net present cost for each size
npc = zeros(1,n);
for i = 1:1:n
    [npc(:,i),costs(:,:,i),mc(:,i)] = design_test(testSystems(SYSINDEX),gen_i_size(i),1,gen_i,TestData,design_day,years);
end
[~,i_min] = min(npc);
if i_min == 1
    gen_i_size = linspace(gen_i_size(1),gen_i_size(2),n);
    npc(1,:) = npc(:,i_min);
    for i = 2:1:n
        [npc(:,i),costs(:,:,i),mc(:,i)] = design_test(testSystems(SYSINDEX),gen_i_size(i),1,gen_i,TestData,design_day,years);
    end
    [~,i_min] = min(npc);
elseif i_min == n
    gen_i_size = linspace(gen_i_size(n-1),gen_i_size(n),n);
    npc(1,:) = npc(:,i_min);
    for i = 1:1:n-1
        [npc(:,i),costs(:,:,i),mc(:,i)] = design_test(testSystems(SYSINDEX),gen_i_size(i),1,gen_i,TestData,design_day,years);
    end
    [~,i_min] = min(npc);
else
    gen_i_size = linspace(gen_i_size(i_min-1),gen_i_size(i_min+1),n);
    npc(1,:) = npc(:,i_min);
    for i = 2:1:n
        [npc(:,i),costs(:,:,i),mc(:,i)] = design_test(testSystems(SYSINDEX),gen_i_size(i),1,gen_i,TestData,design_day,years);
    end
    [~,i_min] = min(npc);
end

% get optimal size
[testSystems(SYSINDEX).Generator,testSystems(SYSINDEX).Costs.Equipment,testSystems(SYSINDEX).Network] = design_resize(testSystems(SYSINDEX).Generator,testSystems(SYSINDEX).Costs.Equipment,testSystems(SYSINDEX).Network,gen_i_size(i_min),1,gen_i);
testSystems(SYSINDEX).Costs.ProjectedMonthlyCosts = mc(:,i_min);
testSystems(SYSINDEX).Costs.NPC = npc(:,i_min);
testSystems(SYSINDEX).Costs.Design = costs(:,:,i_min);
MainScreen1('popupmenu_Callback')
end%ends function optimize_component_size
