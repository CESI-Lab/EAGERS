function coil = cooling_coil_sizing(varargin)
coil = varargin{1};
k = nonzeros((1:length(coil.name))'.*strcmpi(varargin{2},coil.name));
switch coil.type{k}
    case {'DX:TwoSpeed';'DX:SingleSpeed';}%page 805
        if isnan(coil.rated_capacity(k)) %needs sizing
            curve = varargin{3};
            flow = varargin{4};
            T_sup = varargin{5};
            w_sup = varargin{6};
            all_outdoor = varargin{7};
            T_amb = 35;
            if all_outdoor
                T_mix = 26.7;
                w_mix = 0.01115; %T_wb = 19.5, RH = 50.7%
            else
                T_mix = 35;
                w_mix = 0.014; %T_wb = 23.9, RH = 40%
            end
            h_mix = psychometric('Tdb',T_mix,'w',w_mix,'h');
            h_sup = psychometric('Tdb',T_sup,'w',w_sup,'h');

            peak_capacity = flow*(h_mix-h_sup);% kg/s * J/kg  = W
            CapMod = eval_curve(curve,coil.capacity_v_temperature_curve{k},[T_mix,T_amb]); 
            CC_rated = peak_capacity/CapMod;
            if all_outdoor
                CC_rated = max(min(CC_rated,flow/0.00001677),flow/0.00003355);
            else
                CC_rated = max(min(CC_rated,flow/0.00004027),flow/0.00006041);
            end

            rated_slope = (w_mix - w_sup)/(T_mix - T_sup);
            %%find apparatus dewpoint
            T_search = linspace(5,T_sup)';
            W_s = psychometric('Twb',T_search,'w');
            W_coil_linear_extrap = rated_slope*(T_search-T_sup) + w_sup;
            ind = nnz(W_s<W_coil_linear_extrap)+1;
            apparatus_dewpoint = T_search(ind);
            h_ADP = psychometric('Tdb',apparatus_dewpoint,'w',W_s(nnz(W_s<W_coil_linear_extrap)),'h');
            rated_bypass_factor = (h_sup - h_ADP)/(h_mix-h_ADP);
            coil.rated_A0(k,1) = -log(rated_bypass_factor)*flow;

            if isnan(coil.rated_sensible_heat_ratio(k,1))
                h_ADP_search = psychometric('Tdb',T_search_K-273,'w',W_s,'h');
                w_ADP = interp1(h_ADP_search,W_s,h_ADP);
                coil.rated_sensible_heat_ratio(k,1) = min(1,(psychometric('Tdb',T_mix,'w',w_ADP,'h')-h_ADP)/(h_mix - h_ADP));
            end
            switch coil.type{k}
                case 'DX:TwoSpeed'
                    coil.rated_air_flow(k,1) = flow;
                    coil.rated_capacity(k,1) = CC_rated;
                    coil.rated_capacity2(k,1) = CC_rated/2;
                    coil.rated_air_flow2(k,1) = flow/2;
                case 'DX:SingleSpeed'
                    coil.rated_air_flow(k,1) = flow;
                    coil.rated_capacity(k,1) = CC_rated;
            end
        end
    case 'Water:Detailed'
        air_in.T = varargin{3};
        air_out.T = coil.air_outlet_temperature(k);
        air_in.w = varargin{4};
        air_out.w = varargin{4};%%assume same humidity
        air_in.h = psychometric(air_in,'h');
        air_out.h = psychometric(air_out,'h');
        m_a = varargin{5};
        Cp_water = 4186; %J/kg*K
        coil.capacity(k,1) = m_a*(air_in.h-air_out.h);%W
        m_w = coil.capacity(k,1)/(Cp_water*(coil.outlet_water_temperature(k) - coil.inlet_water_temperature(k)));
        coil.tubes_per_row(k,1) = max(3,13750*m_w);
        coil.fin_diameter(k,1) = 0.335*m_a;
        coil.airflow_area(k,1) = 0.44*m_a;
        coil.fin_area(k,1) = 78.5*m_a;
        coil.area_inside(k,1) = 4.4*coil.tube_id(k)*coil.tube_rows(k)*coil.tubes_per_row(k,1);
        coil.area_outside(k,1) = 4.1*coil.tube_od(k)*coil.tube_rows(k)*coil.tubes_per_row(k,1);
        coil.coil_depth(k,1) = coil.tube_spacing(k)*coil.tube_rows(k);
    case 'Water'
        Tw_in0 = coil.water_inlet_temperature(k);
        air_in.T = varargin{3};
        air_out.T = coil.air_outlet_temperature(k);
        air_in.w = varargin{4};
        air_out.w = varargin{4};%%assume same humidity
        air_in.h = psychometric(air_in,'h');
        air_out.h = psychometric(air_out,'h');
        m_a = varargin{5};
        Cp_air = 1006 + 1860*air_in.w; %J/kg*K
        Cp_water = 4186; %J/kg*K
        coil.capacity(k,1) = m_a*(air_in.h-air_out.h);%W
        m_w = coil.capacity(k,1)/(Cp_water*(coil.water_outlet_temperature(k) - coil.water_inlet_temperature(k)));
        %find nominal UA factor
        C_air = Cp_air*m_a;
        C_w = Cp_water*m_w;
        C_min = min(C_air,C_w);
        Z = C_min/max(C_air,C_w);
        max_heat_transfer = C_min*(air_in.T - Tw_in0);
        effectiveness = (air_in.T - air_out.T)*C_air/max_heat_transfer;
        NTU = [0,.001,.01,.1,1,logspace(1,5,5)];
        effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
        for j = 1:1:3
            NTU = linspace(NTU(nnz(effect<effectiveness)),NTU(11-nnz(effect>effectiveness)),10); 
            effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
        end
        NTU = interp1(effect,NTU,effectiveness);
        coil.UA(k,1) = NTU*C_min;
        coil.rated_air_flow(k,1) = m_a;
        coil.rated_water_flow(k,1) = m_w;
        coil.air_inlet_temperature(k,1) = air_in.T;
    otherwise
        disp('need to add missing cooling coil type')
end
end%Ends function cooling_coil_sizing