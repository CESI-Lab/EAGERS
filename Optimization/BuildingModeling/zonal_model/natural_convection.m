function h = natural_convection(normal,dT,nat_buoyancy,method)
switch method
    case 'DOE-2'
%         h = max(.1,1.77*abs(dT).^(1/4));% page 538 https://www.researchgate.net/publication/264308543_COMPARISON_OF_ENERGYPLUS_AND_DOE-2_DETAILED_WINDOW_HEAT_TRANSFER_MODELS
        h = max(.1,1.31*abs(dT).^(1/3));%ASHRAE vertical surface %https://www.nrel.gov/docs/fy12osti/55787.pdf page 18
    case {'Detailed';'BLAST';'TARP';}    
        %% TARP method %A) simple buoyancy
        h = 1.31*abs(dT).^(1/3);%ASHRAE vertical surface
        cos_phi = normal(:,3)./(normal(:,1).^2 + normal(:,2).^2 + normal(:,3).^2).^.5; %portion of surface pointed upwards
        h_tilt_unstable = 9.482*abs(dT).^(1/3)./(7.283-abs(cos_phi));%Walton for tilted surfaces
        h_tilt_stable = 1.81*abs(dT).^(1/3)./(1.382-abs(cos_phi));%Walton for tilted surfaces
        stable  = (normal(:,3)<-1e-3 & dT>0) | (normal(:,3)>1e-3 & dT<0);
        unstable  = (normal(:,3)<-1e-3 & dT<0) | (normal(:,3)>1e-3 & dT>0);
        h(stable) = h_tilt_stable(stable);
        h(unstable) = h_tilt_unstable(unstable);
        h = max(.1,h);
    otherwise
        %% other methods see pg 105 of reference
        disp('need the adaptive convection algorithm')
end
end%ends function natural_convection