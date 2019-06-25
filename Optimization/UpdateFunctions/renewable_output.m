function power = renewable_output(gen,subnet,date,irradiation)
% i is the generator #
%Date is the vector of time points at which you request data
%str chooses between forecast and actual data
n_s = length(date);
n_g = length(gen);
power = zeros(n_s,n_g);
for i = 1:1:n_g
    if strcmp(gen(i).Source,'Renewable') && gen(i).Enabled
        gen_i = gen(i).VariableStruct;
        if isfield(gen_i,'QPform') && isfield(gen_i.QPform,'Electrical')
            location = subnet.Electrical.Location(gen_i.QPform.Electrical.subnetNode);
        else
            location = struct('Latitude',40, 'Longitude',-105, 'TimeZone',-7);
        end
        if strcmp(gen_i.Type,'Solar') 
            [~,~,azimuth,zenith] = solar_calc(location.Longitude,location.Latitude, location.TimeZone,date);
            if strcmp(gen_i.Tracking,'fixed')
                power(:,i) = gen_i.Sizem2 * (irradiation/1000).*cosd(zenith-gen_i.Tilt).*max(0,(cosd(azimuth-gen_i.Azimuth)))*gen_i.Eff;
            elseif strcmp(gen_i.Tracking,'1axis')
                power(:,i) = gen_i.Sizem2*(irradiation/1000).*cosd(zenith-gen_i.Tilt)*gen_i.Eff;
            else
                power(:,i) = gen_i.Sizem2*(irradiation/1000)*gen_i.Eff; %Dual axis
            end
        elseif strcmp(gen_i.Type,'Wind')
            %% need to add wind
        end
    end
end
end%ends function renewable_output