function [air_out,Power] = cooling_coil(building,name,air_in,T_outside,T_treat)
%% Compute the cooling and electric load for DX cooling coils (local AC units not on a central water loop)
coil = building.coils.cooling;
curve = building.curves;
air_density = 1.204; %kg/m^3 %standard air density (1.204 kg/m3) adjusted for the local barometric pressure (standard barometric pressure corrected for altitude, ASHRAE 1997 HOF pg. 6.1).
k = nonzeros((1:length(coil.name))'.*strcmpi(name,coil.name));
air_out = air_in;
spec_heat_supply = 1006 + 1860*air_in.w; %J/kg*K
sensible_load = air_in.m_dot*spec_heat_supply*(air_in.T - T_treat);
if sensible_load>10
    switch coil.type{k}
        case {'DX:SingleSpeed';'DX:TwoSpeed'}%page 805
            c_n = f_index(coil.name{k},building.hvac.components.name);
            if ~isempty(c_n)
                al = building.hvac.components.loop(c_n);%air loop
                bf_rated = bypass_factor_rated(building.hvac.design.Tsupply_c(al),building.hvac.design.w_supply_c(al),building.hvac.design.all_outdoor_flow_cooling(al));
            elseif isfield(building.unitary_heat_cool,'cool_coil')
                uhc_n = f_index(coil.name{k},building.unitary_heat_cool.cool_coil);
                al = building.unitary_heat_cool.loop(uhc_n);
                bf_rated = bypass_factor_rated(building.hvac.design.Tsupply_c(al),building.hvac.design.w_supply_c(al),building.hvac.design.all_outdoor_flow_cooling(al));
            else
                us = f_index(coil.name{k},building.unitary_sys.cool_coil);
                z = building.unitary_sys.zone(us);
                sp = f_index(z,building.setpoints.zone);
                bf_rated = bypass_factor_rated(building.setpoints.cooling_T_des(sp),building.setpoints.cooling_w_des(sp),false);
            end
            coil.rated_A0(k,1) = -log(bf_rated)*coil.rated_air_flow(k)*air_density;
            
            bypass_factor = exp(-coil.rated_A0(k)./air_in.m_dot);
            input(:,1) = psychometric(air_in,'Twb'); 
            input(:,2) = T_outside; %needs to be outdoor wetbulb temperature if it is an evaporative cooled condenser
            
            TotCapTempMod= eval_curve(curve,coil.capacity_v_temperature_curve{k},input);
            EIRTempMod = eval_curve(curve,coil.energy_input_v_temperature_curve{k},input);

            input2(:,1) = air_in.m_dot/(coil.rated_air_flow(k)*air_density); %normalized flow
            TotCapFlowMod = eval_curve(curve,coil.capacity_v_flow_curve{k},input2); 
            EIRFlowMod= eval_curve(curve,coil.energy_input_v_flow_curve{k},input2);
            EIR = (1/coil.rated_COP(k))*EIRTempMod*EIRFlowMod;%energy input ratio
            Q_total = coil.rated_capacity(k)*TotCapTempMod*TotCapFlowMod;
            PLR = min(1,sensible_load/Q_total);
            Q_max = Q_total;
            Speed_ratio = 1;
            if strcmp(coil.type{k},'DX:TwoSpeed')
                TotCapTempMod_2 = eval_curve(curve,coil.capacity_v_temperature_curve2{k},input);
                EIRTempMod_2 = eval_curve(curve,coil.energy_input_v_temperature_curve2{k},input);

                Q_total2 = coil.rated_capacity2(k)*TotCapTempMod_2*TotCapFlowMod;
                EIR_2 = (1/coil.rated_COP2(k))*EIRTempMod_2*EIRFlowMod;%energy input ratio

                Speed_ratio = max(0,(sensible_load - Q_total2)/(Q_total - Q_total2));
                if Speed_ratio<0
                    PLR = min(1,sensible_load/Q_total2);
                    Q_max = Q_total2;
                end
            end
            error = 1;
            while abs(error)>1e-4
                air_out = coil_air_out(air_in,Q_max,bypass_factor,PLR);
                air_out2 = coil_air_out(air_in,Q_max,bypass_factor,PLR+1e-3);
                dT_dPLR = (air_out.T - air_out2.T)/1e-3;
                error = (air_out.T - T_treat);
                if (PLR == 1 && error>0) || (error<0 && PLR == 0)
                    error = 0;
                end
                PLR = max(0,min(1,PLR + error/dT_dPLR));
            end
            if strcmp(coil.type{k},'DX:TwoSpeed') && Speed_ratio>0
                revised_sensible_load = PLR*Q_total;
                Speed_ratio = max(0,(revised_sensible_load - Q_total2)/(Q_total - Q_total2));
                Power = Speed_ratio*Q_total*EIR + (1-Speed_ratio)*Q_total2*EIR_2;
                heat_rejected = Speed_ratio*Q_total*(1+EIR) + (1-Speed_ratio)*Q_total2*(1+EIR_2);
            else
                PartLoadFrac = max(.7,eval_curve(curve,coil.part_load_curve{k},PLR));
                RTF = PLR./PartLoadFrac;%runtime fraction
                if strcmp(coil.type{k},'DX:TwoSpeed')
                    Power = Q_total2*EIR_2*RTF;
                    heat_rejected = Q_total2*(1+EIR_2);
                else
                    Power = Q_total*EIR*RTF;
                    heat_rejected = Q_total*(1+EIR);
                end
            end
    end
else
    Power = 0;
end
end%Ends function cooling_coil

function air_out = coil_air_out(air_in,Q_total,bypass_factor,PLR)
%% Determin exit air conditions based on operating status of the coil

%find apparatus dewpoint
h_ADP = air_in.h - (Q_total/air_in.m_dot)/(1-bypass_factor);% J/kg
T_search = linspace(5,40,20)';
w_ADP = interp1(psychometric('Twb',T_search,'h'),psychometric('Twb',T_search,'w'),h_ADP);

air_out.m_dot = air_in.m_dot;
air_out.h = air_in.h - Q_total*PLR/air_out.m_dot;%J/kg
if air_in.w<w_ADP %dry coil
    sensible_heat_ratio = 1;
    air_out.w = air_in.w;
else
    sensible_heat_ratio = min(1,(psychometric('Tdb',air_in.T,'w',w_ADP,'h')-h_ADP)/(air_in.h - h_ADP));
    h_Tin_w_out = air_in.h - (1-sensible_heat_ratio)*(air_in.h - air_out.h);
    air_out.w = psychometric('Tdb',air_in.T,'h',h_Tin_w_out,'w'); 
end
spec_heat = 1006 + 1860*(air_in.w+air_out.w)/2; %J/kg*K
air_out.T = air_in.T - sensible_heat_ratio*Q_total*PLR/(air_out.m_dot*spec_heat);
end%Ends function coil_air_out


function bf = bypass_factor_rated(T_supply,w_supply,all_outdoor)
%% determine bypass factor based on rated dew point
if all_outdoor
    T_mix = 26.7;
    w_mix = 0.01115; %T_wb = 19.5, RH = 50.7%
else
    T_mix = 35;
    w_mix = 0.014; %T_wb = 23.9, RH = 40%
end
h_mix = psychometric('Tdb',T_mix,'w',w_mix,'h');
h_sup = psychometric('Tdb',T_supply,'w',w_supply,'h');

rated_slope = (w_mix - w_supply)/(T_mix - T_supply);
%%find apparatus dewpoint
T_search = linspace(5,T_supply,20)';
W_s = psychometric('Twb',T_search,'w');
W_coil_linear_extrap = rated_slope*(T_search-T_supply) + w_supply;
ind = min(20,nnz(W_s<W_coil_linear_extrap)+1);
apparatus_dewpoint = T_search(ind);

h_ADP = psychometric('Tdb',apparatus_dewpoint,'w',W_s(ind),'h');
bf = (h_sup - h_ADP)/(h_mix-h_ADP);
end%Ends function bypass factor