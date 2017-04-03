% RUN SCRIPT FROM genie-matlab folder!
% OR load load('Seiter_TOC_Obs_03042017.mat')
load('TOC_Seiteretal2004.asc')
% change -9999 to NaN
TOC_Seiteretal2004(TOC_Seiteretal2004==-9999)=NaN;
% ALSO adjust outlier 30100000 in 134, 285 to 0.3
[r,c] = size(TOC_Seiteretal2004)
test = NaN(2,c);
TOC_NEW = [test;TOC_Seiteretal2004];
TOC_NEW = [TOC_NEW, test];
TOC_NEW = [TOC_NEW; test];
TOC_NEW_tr = TOC_NEW.'
[dum_lon, dum_lat] = get_grid_genie_lonlatedges(36,36,0);
[dum_lon_orig, dum_lat_orig] = get_grid_genie_lonlatedges(361,181,0);
[zo fao] = make_regrid_2d(dum_lon_orig,dum_lat_orig,TOC_NEW_tr,dum_lon,dum_lat,true)