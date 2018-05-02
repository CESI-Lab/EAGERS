function design_grid_search
global testSystems SYSINDEX GENINDEX TestData

project = testSystems(SYSINDEX);
gen_i = testSystems(SYSINDEX).Generator(GENINDEX);
gen_i_size = [(0.2*gen_i.Size : 0.2*gen_i.Size : gen_i.Size),1.5*gen_i.Size,2*gen_i.Size];
handles = guihandles;
if get(handles.DesignDay,'Value') == 1
    design_day = true;
else
    design_day = false;
end
years = str2double(get(handles.NPC_Years,'String'));%choose years in GUI
% find net present cost for each size
npc = zeros(size(gen_i_size));
for i = 1:1:length(gen_i_size)
    [npc(:,i),costs(:,:,i),mc(:,i)] = design_test(project,gen_i_size(i),GENINDEX,TestData,design_day,years);
end

% get optimal size
[~,index_minimum] = min(npc);
[testSystems(SYSINDEX).Generator,testSystems(SYSINDEX).Costs.Equipment] = design_resize(testSystems(SYSINDEX).Generator,testSystems(SYSINDEX).Costs.Equipment,gen_i_size(index_minimum),GENINDEX);
testSystems(SYSINDEX).Costs.ProjectedMonthlyCosts = mc(:,index_minimum);
testSystems(SYSINDEX).Costs.NPC = npc(:,index_minimum);
testSystems(SYSINDEX).Costs.Design = costs(:,:,index_minimum);
% [testSystems(SYSINDEX).Costs.Financial.irr,testSystems(SYSINDEX).Costs.Financial.payback] = FinancialMetric(testSystems(SYSINDEX).Costs.ProjectedMonthlyCosts,testSystems(SYSINDEX).Costs.ProjectedMonthlyCosts,(1+test_sys(i_ts).Costs.DiscountRate/100));
MainScreen1('popupmenu_Callback')