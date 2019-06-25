function [zip_code,state] = get_zip_code(lat,lon)
global Model_dir
load(fullfile(Model_dir,'Data','ZipMap.mat'))
[~,index] = min(abs(ZipInfo(:,2) - lat) + abs(ZipInfo(:,3) - lon));
zip_code = ZipInfo(index,1);
load(fullfile(Model_dir,'Data','ZipState.mat'))
[~,index] = min(abs(ZipCodeInfo.ZipCode - zip_code));

state_name = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
         'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
         'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
         'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
state_abbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};            
abbrev = ZipCodeInfo.State(index);
state = state_name{nonzeros((1:length(state_abbrev))'.*strcmp(abbrev,state_abbrev))};
end%Ends function get_zip_code
