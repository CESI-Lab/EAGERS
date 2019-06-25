function coil = heating_coil_sizing(coil,name,T_air,w_air,m_air)
%used for water coils only
k = nonzeros((1:length(coil.name))'.*strcmpi(name,coil.name));
Tw_in0 = coil.water_inlet_temperature(k);
air_in.T = T_air;
air_out.T = coil.air_outlet_temperature(k);
air_in.w = w_air;
air_out.w = w_air;%%assume same humidity
air_in.h = psychometric(air_in,'h');
air_out.h = psychometric(air_out,'h');
m_a = m_air;
Cp_air = 1006 + 1860*air_in.w; %J/kg*K
Cp_water = 4186; %J/kg*K
coil.capacity(k,1) = m_a*(air_out.h-air_in.h);%W
m_w = coil.capacity(k,1)/(Cp_water*(coil.water_inlet_temperature(k) - coil.water_outlet_temperature(k)));

C_air = Cp_air*m_a;
C_w = Cp_water*m_w;
C_min = min(C_air,C_w);
Z = C_min/max(C_air,C_w);
max_heat_transfer = C_min*(Tw_in0 - air_in.T);
effectiveness = (air_out.T - air_in.T)*C_air/max_heat_transfer;
NTU = [0,.001,.01,.1,1,logspace(1,5,5)];
effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
for j = 1:1:3
    nn = nnz(effect<effectiveness);
    NTU = linspace(NTU(nn),NTU(11-nn),10); 
    effect = 1 - exp((exp(-Z*NTU.^.78)-1)./(Z*NTU.^(-.22)));
end
NTU = interp1(effect,NTU,effectiveness);
coil.UA(k,1) = NTU*C_min;
coil.rated_air_flow(k,1) = m_a;
coil.rated_water_flow(k,1) = m_w;
if isnan(coil.air_inlet_temperature(k,1))
    coil.air_inlet_temperature(k,1) = air_in.T;
end
end%Ends function heating_coil_sizing