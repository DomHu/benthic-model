% Re-gridding modern TOC wt% observations from Seitert et al. 2004
% either in TOC_Seiteretal2004.asc where I need to do the adjustments 
% which are commented out below 
% or load Seiter_TOC_Obs_03042017.mat where this has been done already

% run within genie-matlab folder!

load('Seiter_TOC_Obs_03042017.mat');
% % % % change -9999 to NaN
% % % TOC_Seiteretal2004(TOC_Seiteretal2004==-9999)=NaN;
% % % % ALSO adjust outlier 30100000 in 134, 285 to 0.3
% % % [r,c] = size(TOC_Seiteretal2004)
% % % test = NaN(2,c);
% % % TOC_NEW = [test;TOC_Seiteretal2004];
% % % TOC_NEW = [TOC_NEW, test];
% % % TOC_NEW = [TOC_NEW; test];
% % % 
% one column too much: Avg first and last column - omit NaN
help=horzcat(TOC_NEW(:,1), TOC_NEW(:,end));
Nmean=nanmean(help,2);
% use this as first column
TOC_NEW(:,1)=Nmean;
% delete the last column
TOC_NEW = TOC_NEW(:,1:end-1);
TOC_NEW_tr = TOC_NEW.';

[dum_lon, dum_lat] = get_grid_genie_lonlatedges(36,36,-180);

%[dum_lon_orig, dum_lat_orig] = get_grid_genie_lonlatedges(361,181,0);
% create edges of observations manually:
long_obs=(-180.5:1:179.5);
lat_obs=(-90.5:1:90.5);
lat_obs(1)=-90;
lat_obs(end)=90;

[zo fao] = make_regrid_2d(long_obs,lat_obs,TOC_NEW_tr,dum_lon,dum_lat,true);

dlmwrite('zo_Seiter_offset1905.txt',zo.')
% zo or go_z is what I need!