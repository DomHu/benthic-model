function DOM_Seiter_make_lat_lon_obs()
Seiter_data=xlsread('Seiter_on_GENIE_72x72.xlsx','Sheet1','C154:BV225');
Seiter_lat=xlsread('Seiter_on_GENIE_72x72.xlsx','Sheet1','B154:B225');
Seiter_lon=xlsread('Seiter_on_GENIE_72x72.xlsx','Sheet1','C153:BV153')';

% just 70 rows because in the very south just NaN
sz = size(Seiter_data);
rep_lon= kron(Seiter_lon, ones(sz(1),1));
rep_lat=kron(ones(sz(2),1), Seiter_lat);

rep_lon_lat = cat(2, rep_lon, rep_lat);
Seiter_overlay_data= cat(2, rep_lon_lat, Seiter_data(:));

end