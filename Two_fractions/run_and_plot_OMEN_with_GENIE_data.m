function [STATM] = run_and_plot_OMEN_with_GENIE_data(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME, POUTVAR, PKPARAM, PK2SCALING)
%
%   *******************************************************************   %
%   *** sedgem 2-D (LON-LAT) DATA PLOTTING ****************************   %
%   *******************************************************************   %
%
%   run_and_plot_OMEN_with_GENIE_data(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME, POUTVAR)
%   plots the SEDGEM 2-D netCDF data file 'fields_sdegem_2d.nc' and takes 15 arguments:
%
%   PEXP1 [STRING] (e.g. 'preindustrial_spinup')
%   --> the (first) experiment name
%   PEXP2 [STRING] [OPTIONAL] (e.g. 'enhanced_export')
%   --> the experiment name of 2nd, optional, netCDF file
%   --> leave EXP2 blank, i.e., '', for no second netCDF file
%   PVAR1 [STRING] (e.g. 'zox')
%   -->  the name of the variable to be plotted
%   --> also part of the title and file-name
%   PVAR2 [STRING] [OPTIONAL] (e.g. 'sed_opal')
%   --> id the name of 2nd, optional, variable
%   PT1 [REAL] [OPTIONAL] (e.g. 0.0)
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PT2 [REAL] [OPTIONAL] (e.g. 0.0)
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PIK [INTEGER] [OPTIONAL] (e.g. 0)
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PMASK [STRING] [OPTIONAL] (e.g. '')
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PCSCALE [REAL] (e.g. 1.0)
%   --> the scale factor for the plot
%       e.g., to plot as weight fraction (instead of wt%), enter: 0.01
%             to plot negative values, enter: -1
%   --> the plot is auto-scaled if a value of zero (0.0) is entered
%   PCMIN [REAL] (e.g. 0.0)
%   --> the minimum scale value
%   PCMAX [REAL] (e.g. 100.0)
%   --> the maxumum scale value
%   PCN [INTEGER] (e.g. 20)
%   --> the number of (contor) intervals between minimum and maximum
%       scale values
%   PDATA [STRING] (e.g. 'obs_d13C.dat')
%   --> the filename containing any overlay data set,
%       which must be formatted as (space) seperated columns of:
%       lon, lat, value, label
%       or, of the option (below) data_ij = 'y', then as:
%       i, j, value, label
%   --> the full filename must be give, including any extensions
%   --> leave PDATA blank, i.e., '', for no overlay data
%   POPT [STRING] (e.g., 'plotting_config_2')
%   --> the string for an alternative plotting parameter set
%   --> if an empty (i.e., '') value is passed to this parameter
%       then the default parameter set is used
%   PNAME [STRING] (e.g., 'my_plot')
%   --> the string for an alternative filename
%   --> if an empty (i.e., '') value is passed to this parameter
%       then a filename is automatically generated
%
%   POUTVAR [INTEGER] (e.g. 1, 2, ...)
%   --> number of the output to be calculated:
%   -->  1: zox
%   -->  2: zSO4
%   -->  3: depth integrated total OM oxidation rate (Cox total) [mol cm-2 yr-1]
%   -->  4: depth integrated aerobic OM oxidation rate (Cox aerobic) [mol cm-2 yr-1]
%   -->  5: depth integrated anaerobic (SO4 reduction at the moment) OM oxidation rate (Cox SO4) [mol cm-2 yr-1]
%   -->  6: fraction of O2/total depth integrated OM oxidation
%   -->  7: fraction of SO4/total depth integrated OM oxidation
%   -->  8: total TOC wt% preserved at a specific depth (e.g. zinf = 100cm - to be specified in benthic_test.OMEN_with_GENIE_input)
%   -->  9: labile TOC1 wt% preserved at a specific depth (e.g. zinf = 100cm)
%   --> 10: refractory TOC2 wt% preserved at a specific depth (e.g. zinf = 100cm)
%   --> 11: SWI-flux O2
%   --> 12: SWI-flux SO4
%   --> 13: SWI-flux H2S
%   --> 14: SWI-flux PO4
%   --> 15: SWI-flux DIC
%   --> 16: SWI-flux ALK       
%   --> 17: mean OM wt% in upper x cm
%   PKPARAM [STRING] (e.g., 'boudreau1997')
%   --> the string for the k parameterisation scheme
%   -->  options are 
%   -->  'boudreau1997'
%   -->  'tromp1995'
%   -->  'stolpovsky2016'
%   -->  'boudreau1997fPOC'
%   -->  'invariant' here, global values must be specified in
%         benthic_test.OMEN_with_GENIE_input(...)
%   --> PK2SCALING [REAL] scaling factor for calculating k2 = PK2SCALING * k1
%
%   Example
%       1)  run_and_plot_OMEN_with_GENIE_data('0804_13_worjh2_closed_No_PO4_Full_OMEN_k_Tromp','','oxygen penetration depth','',0.0,0.0,0,'',1.0,0.0,10.0,20,'','','',1)
%           will plot the oxygen penetration depth (zox),
%           between 0 and 10 cm in 20 contour intervals
%       2) run_and_plot_OMEN_with_GENIE_data('1404_02_worjh2_closed_with_PO4_Full_OMEN_zbio_k_0.01','','frac_of_aerobic_Cox','',0.0,0.0,0,'',1.0,0.0,80.0,10,'','','',6)
%       3) run_and_plot_OMEN_with_GENIE_data('1404_02_worjh2_closed_with_PO4_Full_OMEN_zbio_k_0.01','','Total_Cox_rate','',0.0,0.0,0,'',1e-6,0.0,400.0,10,'','','',3)
%       4) run_and_plot_OMEN_with_GENIE_data('1404_02_worjh2_closed_with_PO4_Full_OMEN_zbio_k_0.01','','TOC wtpc at 100cm','',0.0,0.0,0,'',1.0,0.0,5.0,10,'','','',8)
%       5) run_and_plot_OMEN_with_GENIE_data(exp_1,'','Mean OM wt% in upper 10cm','',0.0,0.0,0,'',1.0,0.0,5.0,10,'','','',17, k_parametr)
%   *******************************************************************   %

% *********************************************************************** %
% ***** HISTORY ********************************************************* %
% *********************************************************************** %
%
%   17/04/06: used plot_fields_sedgem_2d as starting point
%
% *********************************************************************** %
%%

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% close open windows
% NOTE: don't clear variable space here ...
close all;

% set paths
my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
addpath([ my_dir '/../Two_fractions'])

% load plotting options
if isempty(POPT), POPT='plot_fields_settings'; end
eval(POPT);
% set axes
lat_min = -090;
lat_max = +090;
lon_min = plot_lon_origin;
lon_max = lon_min+360;
lon_offset = 0;
% set passed parameters
exp_1 = PEXP1;
exp_2 = PEXP2;
dataid_1 = PVAR1;
dataid_2 = PVAR2;
dummy_1 = PT1;
dummy_2 = PT2;
dummy_3 = PIK;
maskid = PMASK;
data_scale = PCSCALE;
con_min = PCMIN;
con_max = PCMAX;
con_n = PCN;
overlaydataid = PDATA;
altfilename = PNAME;
plotvar = POUTVAR;
k_parametr = PKPARAM;
k2scaling = PK2SCALING;
% DOM: data needed from GENIE:
Nitrogen=false;
data_fsed_POC='fsed_POC      ';
%dataid_1=strtrim(data_fsed_POC);
data_ocn_O2='ocn_O2        ';
data_ocn_NO3='ocn_NO3       ';
data_ocn_NH4='ocn_NH4       ';
data_ocn_SO4='ocn_SO4       ';
data_ocn_H2S='ocn_H2S       ';
data_ocn_PO4='ocn_PO4       ';
data_ocn_ALK='ocn_ALK       ';
data_ocn_DIC='ocn_DIC       ';
data_burial_rate_OMEN='misc_OMEN_bur ';
% now use actual burial flux from previous time-step (also saved as netcdf now)
% substitute last data_GENIE from below (data_fburial_POC='fsed_POC) with
% this value

% This is sedimentation flux
data_fburial_CaCO3='fsed_CaCO3    ';
data_fburial_det='fsed_det      ';
data_fburial_ash='fsed_ash      ';
% rather use the preservation flux: (units = "mol cm-2 yr-1")
%data_fburial_CaCO3='fburial_CaCO3 ';
% data_fburial_det='fburial_det   ';
data_fburial_POC='fsed_POC      ';
% data_fburial_ash='fburial_ash   ';

data_ben_temp='ocn_temp      ';
data_depth='grid_topo     ';
data_is_frac2='fsed_POC_frac2';

if(Nitrogen)
    data_GENIE = [data_fsed_POC; data_ocn_O2; data_ocn_NO3; data_ocn_NH4; data_ocn_SO4; data_ocn_H2S; data_ocn_PO4; data_ocn_DIC; data_ocn_ALK; data_fburial_CaCO3; data_fburial_det; data_fburial_ash; data_ben_temp; data_depth; data_is_frac2; data_burial_rate_OMEN];
%    data_GENIE = [data_fsed_POC; data_ocn_O2; data_ocn_NO3; data_ocn_NH4; data_ocn_SO4; data_ocn_H2S; data_ocn_PO4; data_ocn_DIC; data_ocn_ALK; data_fburial_CaCO3; data_fburial_det; data_fburial_ash; data_ben_temp; data_depth; data_is_frac2; data_burial_rate_OMEN];
else
    data_GENIE =[data_fsed_POC; data_ocn_O2; data_ocn_SO4; data_ocn_H2S; data_ocn_PO4; data_ocn_DIC; data_ocn_ALK; data_fburial_CaCO3; data_fburial_det; data_fburial_ash; data_ben_temp; data_depth; data_is_frac2; data_burial_rate_OMEN];
%    data_GENIE =[data_fsed_POC; data_ocn_O2; data_ocn_SO4; data_ocn_H2S; data_ocn_PO4; data_ocn_DIC; data_ocn_ALK; data_fburial_CaCO3; data_fburial_det; data_fburial_ash; data_ben_temp; data_depth; data_is_frac2; data_burial_rate_OMEN];
end
n_data_GENIE = size(data_GENIE,1);

% set default data scaling
if data_scale == 0.0
    data_scale = 1.0;
    clear con_min;
    clear con_max;
    con_n = 10;
end
% set global flag if no alt plotting scale is set
% NOTE: catch possibility of one axis being set, but the other @ default
%       (min and max with indetical values)
if ((plot_lat_min == plot_lat_max) && (plot_lon_min == plot_lon_max)),
    plot_global = true;
    plot_xy_scaling = 1.0;
else
    plot_global = false;
    if (plot_lat_min == plot_lat_max),
        plot_lat_min = lat_min;
        plot_lat_max = lat_max;
    end
    if (plot_lon_min == plot_lon_max),
        plot_lon_min = lon_min;
        plot_lon_max = lon_max;
    end
    plot_xy_scaling = ((plot_lat_max - plot_lat_min)/(lat_max - lat_min)) / ((plot_lon_max - plot_lon_min)/(lon_max - lon_min));
end
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% set function name
str_function = 'plot-fields-sedgem-2d';
% plot format
if ~isempty(plot_format), plot_format_old='n'; end
% plotting paths
addpath(library_path1);
if (plot_format_old == 'n'),
    addpath(library_path2);
    addpath(library_path3);
end
%
% *** DEFINE COLORS ***************************************************** %
%
% define contonental color
color_g = [0.75 0.75 0.75];
% define no-data color
color_b = [0.00 0.00 0.00];
%
% *********************************************************************** %

% *********************************************************************** %
% *** INITIALIZE ARRAYS ************************************************* %
% *********************************************************************** %
%
xm = [];
ym = [];
zm = [];
lonm = [];
lone = [];
lonw = [];
latm = [];
latn = [];
lats = [];
laym = [];
layt = [];
layb = [];
rawdata=[];
data_1=[];
data_2=[];
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPEN netCDF DATA FILE ********************************************* %
% *********************************************************************** %
%
% open netCDF file -- test for 'experiment format' or not
% NOTE: other format is indicated by '.nc' extension as experiment ID
if strcmp(exp_1(end-2:end),'.nc'),
    ncid_1=netcdf.open(exp_1,'nowrite');
else
    ncid_1=netcdf.open([data_path '/' exp_1 '/sedgem/fields_sedgem_2d.nc'],'nowrite');
end
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_1);
% load mask data
% NOTE: flip in j-direction to make consistent with netCDF grid
maskfile = maskid;
if ~isempty(maskid)
    mask = load(maskfile,'-ascii');
    mask = flipdim(mask,1);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% load grid data
varid  = netcdf.inqVarID(ncid_1,'grid_topo');
grid_topo(:,:) = netcdf.getVar(ncid_1,varid);
% flip array around diagonal to give (j,i) array orientation
grid_topo = grid_topo';
% calculate grid dimensions
varid  = netcdf.inqDimID(ncid_1,'lat');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
jmax = dimlen;
varid  = netcdf.inqDimID(ncid_1,'lon');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
imax = dimlen;
% load remaining grid information
varid  = netcdf.inqVarID(ncid_1,'lat');
grid_lat = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon');
grid_lon = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonm latm] = meshgrid(grid_lon,grid_lat);
varid  = netcdf.inqVarID(ncid_1,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonw lats] = meshgrid(grid_lon_edges(1:imax),grid_lat_edges(1:jmax));
[lone latn] = meshgrid(grid_lon_edges(2:imax+1),grid_lat_edges(2:jmax+1));
% Non-uniform lat grid
if (plot_equallat == 'n'),
    lat_max = sin(pi*lat_max/180.0);
    lat_min = sin(pi*lat_min/180.0);
    latn = sin(pi*latn/180.0);
    lats = sin(pi*lats/180.0);
    plot_lat_max = sin(pi*plot_lat_max/180.0);
    plot_lat_min = sin(pi*plot_lat_min/180.0);
end
%initialize dummy topography
topo = zeros(jmax,imax);
layb = zeros(jmax,imax);
%i=1 longitude grid origin
grid_lon_origin = grid_lon_edges(1);
% internal mask
if ~isempty(plot_mask_netcdf)
    varid  = netcdf.inqVarID(ncid_1,plot_mask_netcdf);
    rawdata = netcdf.getVar(ncid_1,varid);
    grid_mask(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD DIFFERENT DATASETS         *********************************** %
% *** was SET PRIMARY GRIDDED DATASET *********************************** %
% *********************************************************************** %
%
% check that the variable name (dataid_1) exists - sets varid to n, which is then used to load data 
% Dom: Check for the diffirent data to load and then load all of them
for(k = 1:n_data_GENIE)
    GENIE_dataid = strtrim(data_GENIE(k,:));
    varid = [];
    while isempty(varid)
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
            if strcmp(varname,GENIE_dataid)
                varid = n;
            end
        end
        if isempty(varid)
            disp('   > WARNING: Variable name must be one of the following;');
            for n = 0:nvars-1,
                [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
                varname
            end
            dataid = input('   > Variable name: ','s');
        end
    end
    % load data
    % flip array around diagonal to give (j,i) array orientation
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
    rawdata(:,:,k) = netcdf.getVar(ncid_1,varid);
    rawdata(1:imax,1:jmax,k) = rawdata(1:imax,1:jmax,k);
    data_1(1:jmax,1:imax,k) = rawdata(1:imax,1:jmax,k)';
end 
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET ALTERNATIVE GRIDDED DATASET *********************************** %
% *********************************************************************** %
%
% *** ALT EXPERIMENT **************************************************** %
%
% open new netCDF file if necessary, else reuse 1st netCDF dataset
if ~isempty(exp_2)
    % open netCDF file -- test for 'experiment format' or not
    % NOTE: other format is indicated by '.nc' extension as experiment ID
    if strcmp(exp_2(end-2:end),'.nc'),
        ncid_2=netcdf.open(exp_2,'nowrite');
    else
        ncid_2=netcdf.open([data_path '/' exp_2 '/sedgem/fields_sedgem_2d.nc.nc'],'nowrite');
    end
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_2);
else
    ncid_2 = ncid_1;
end
%
% *** ALT DATA FIELD **************************************************** %
%
if ~isempty(dataid_2)
    % check that the variable name exists
    varid = [];
    while isempty(varid)
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
            if strcmp(varname,dataid_2)
                varid = n;
            end
        end
        if isempty(varid)
            disp('   > WARNING: Variable name must be one of the following;');
            for n = 0:nvars-1,
                [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
                varname
            end
            dataid = input('   > Variable name: ','s');
        end
    end
end
%
% *** SET DATA ********************************************************** %
%
if (~isempty(exp_2)) || (~isempty(dataid_2))
    % load data
    % flip array around diagonal to give (j,i) array orientation
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    rawdata = netcdf.getVar(ncid_2,varid);
    rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
    data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    data_anomoly = 'y';
else
    data_2(:,:,:) = 0.0;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET OUTPUT FILESTRING ********************************************* %
% *********************************************************************** %
%
% create an output filestring for data and plot saving
if ~isempty(exp_2),
    filename = [exp_1, '_MINUS_' exp_2, '.', dataid_1];
    if ~isempty(dataid_2),
        filename = [exp_1, '_MINUS_' exp_2, '.', dataid_1, '_MINUS_', dataid_2];
    end
elseif ~isempty(dataid_2)
    filename = [exp_1, '.', dataid_1, '_MINUS_', dataid_2];
else
    filename = [exp_1, '.', dataid_1];
end
if ~isempty(overlaydataid),
    filename = [filename, '_VS_', overlaydataid];
    if (data_anomoly == 'y'),
        filename = [filename, '.ANOM'];
    end
    if (data_only == 'y'),
        filename = [filename, '.DO'];
    end
end
if ~isempty(maskid),
    filename = [filename, '.', maskid];
end
if (~isempty(altfilename)), filename = altfilename; end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FILTER & PROCESS RAW DATA ***************************************** %
% *********************************************************************** %
%
% set ocean grid value to give white when plotted
if strcmp(dataid_1,'grid_mask')
    data_1 = NaN;
end
if strcmp(dataid_2,'grid_mask')
    data_2 = NaN;
end
%
xm = lonm;
ym = latm;
data = data_1 - data_2;
data = data - data_offset;
zm = data;
%OMEN_results = zeros(jmax,imax,k);
% filter gridded data: Dom calculate values in zm to be plotted further
% down - this should be the results from OMEN (so call OMEN before)
n = 0;
for i = 1:imax,
    for j = 1:jmax,
        if grid_topo(j,i) >= 0
            zm(j,i,:) = NaN;
%            OMEN_results(j,i,:) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        elseif (zm(j,i) < -1.0E6) || (zm(j,i) > 0.9E30) || isnan(zm(j,i))
            zm(j,i,:) = NaN;
%            OMEN_results(j,i,:) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
            % if mask selected and cell is masked out => mark to plot in white
            if ~isempty(plot_mask_netcdf)
                if grid_mask(j,i) == 0.0
                    topo(j,i) = -1.0;
                    layb(j,i) = +1.0;
                end
            end
        else
            if ~isempty(maskid)
                if mask(j,i) == 0
                    zm(j,i,:) = NaN;
%                    OMEN_results(j,i,:) = NaN;
                end
            end
            if data_log10 == 'y'
                if (zm(j,i) > 0.0)
                    zm(j,i,:) = log10(zm(j,i)/data_scale);
                else
                    zm(j,i,:) = NaN;
%                    OMEN_results(j,i,:) = NaN;
               end
            else    % Dom: here scaling is done - also call here OMEN, as these are the ocean grid cells!
                bc=zm(j,i,:);
                i
                j
                test=benthic_test.OMEN_with_GENIE_input(bc, Nitrogen, k_parametr, k2scaling);
%   -->  1: zox
%   -->  2: zSO4
%   -->  3: depth integrated total OM oxidation rate (Cox total)
%   -->  4: depth integrated aerobic OM oxidation rate (Cox aerobic)
%   -->  5: depth integrated anaerobic (SO4 reduction at the moment) OM oxidation rate (Cox SO4)
%   -->  6: fraction of O2/total depth integrated OM oxidation
%   -->  7: fraction of SO4/total depth integrated OM oxidation
%   -->  8: total TOC wt% preserved at zinf (100cm)
%   -->  9: labile TOC1 wt% preserved at zinf (100cm)
%   --> 10: refractory TOC2 wt% preserved at zinf (100cm)
%   --> 11: SWI-flux O2
%   --> 12: SWI-flux SO4
%   --> 13: SWI-flux H2S
%   --> 14: SWI-flux PO4
%   --> 15: SWI-flux DIC
%   --> 16: SWI-flux ALK
%   --> 17: mean OM concentration in upper x cm

                switch plotvar
                case 1
                    OMEN_result(j,i) = test.zox;
                case 2
                    OMEN_result(j,i) = test.zso4;
                case 3
                    OMEN_result(j,i) = test.Cox_rate_total;
                case 4
                    OMEN_result(j,i) = test.Cox_rate_aerobic;
                case 5
                    OMEN_result(j,i) = test.Cox_rate_sulfred;
                case 6
                    OMEN_result(j,i) = test.Cox_perc_aerobic;
                case 7
                    OMEN_result(j,i) = test.Cox_perc_sulfred;
                case 8
                    OMEN_result(j,i) = test.C_zinf_wtpc;
                case 9
                    OMEN_result(j,i) = test.C1_zinf_wtpc;
                case 10
                    OMEN_result(j,i) = test.C2_zinf_wtpc;
                case 11
                    OMEN_result(j,i) = test.flxswiO2;
                case 12
                    OMEN_result(j,i) = test.flxswiSO4;
                case 13
                    OMEN_result(j,i) = test.flxswiH2S;
                case 14
                    OMEN_result(j,i) = test.flxswi_P;
                case 15
                    OMEN_result(j,i) = test.flxswiDIC;
                case 16
                    OMEN_result(j,i) = test.flxswiALK;
                 case 17
                    OMEN_result(j,i) = test.Mean_OM;
                case 18
                    OMEN_result(j,i) = test.Cox_rate_total_xcm;
                 case 19
                    OMEN_result(j,i) = test.Cox_perc_aerobic_xcm;                   
               	otherwise
                    error('Error. Unknown OMEN output specified.')

                end
                
%                OMEN_results(j,i,1) = OMEN_zox;
                zm(j,i) = OMEN_result(j,i)/data_scale;
%                zm(j,i,:) = zm(j,i)/data_scale;
            end
            if ~isnan(zm(j,i)), n = n + 1; end
            topo(j,i) = -1.0;
            layb(j,i) = +1.0;
        end
    end
end
nmax = n;

%%%%% plot zox vs depth
if (plotvar==1)
    zox = OMEN_result(:)';
    wdepth=zm(:,:,end)';
    wdepth=wdepth(:)';
   scatter(zox, wdepth) 
end

%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD (OPTIONAL) OVERLAY DATA ************************************** %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % load overlay datafile
    overlaydatafile = [overlaydataid];
    fid = fopen(overlaydatafile);
    if (data_shapecol == 'n'),
        % lon, lat, value, LABEL
        C = textscan(fid, '%f %f %f %s', 'CommentStyle', '%');
        overlaydata_raw = cell2mat(C(1:3));
        CC = C(4);
        overlaylabel_raw = char(CC{1}(:));
    else
        % lon, lat, value, LABEL, SHAPE, EDGE COLOR, FILL COLOR
        C = textscan(fid, '%f %f %f %s %s %s %s', 'CommentStyle', '%');
        overlaydata_raw = cell2mat(C(1:3));
        CC = C(4);
        overlaylabel_raw = char(CC{1}(:));
        CC = C(5);
        overlaydata_shape = char(CC{1}(:));
        CC = C(6);
        overlaydata_ecol = char(CC{1}(:));
        CC = C(7);
        overlaydata_fcol = char(CC{1}(:));
    end
    fclose(fid);
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    % create (i,j) from (lon,lat) and vice versa (depending on data input type)
    if (data_ijk == 'n')
        % convert (lon,lat) overlay data to (i,j)
        % NOTE: function 'calc_find_ij' takes input in order: (lon,lat)
        %       i.e., the same as the raw overlay data, which is (lon,lat) (i.e., (i,j)) format
        % NOTE: !!! gridded data is (j,i) !!!
        overlaydata_ij(:,:) = zeros(size(overlaydata_raw(:,:)));
        for n = 1:nmax,
            overlaydata_ij(n,1:2) = calc_find_ij(overlaydata_raw(n,1),overlaydata_raw(n,2),grid_lon_origin,imax,jmax);
        end
        overlaydata_ij(:,3) = overlaydata_raw(:,3);
    else
        % convert (i,j) overlay data to (lon,lat)
        % NOTE: save (i,j) data first
        overlaydata_ij(:,:) = overlaydata_raw(:,:);
        overlaydata_raw(:,1) = grid_lon_origin + 360.0*(overlaydata_raw(:,1) - 0.5)/jmax;
        overlaydata_raw(:,2) = 180.0*asin(2.0*(overlaydata_raw(:,2) - 0.5)/jmax - 1.0)/pi;
    end
    % remove data in land cells
    if (data_land == 'n')
        for n = 1:nmax,
            if isnan(zm(overlaydata_ij(n,2),overlaydata_ij(n,1)))
                overlaydata_raw(n,3) = NaN;
                overlaydata_ij(n,3)  = NaN;
            end
        end
        overlaylabel_raw(isnan(overlaydata_raw(:,3)),:) = [];
        overlaydata_raw(isnan(overlaydata_raw(:,3)),:) = [];
        overlaydata_ij(isnan(overlaydata_ij(:,3)),:) = [];
    end
    % update value of nmax
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    %
    overlaylabel(:,:) = overlaylabel_raw(:,:);
    % convert lat to sin(lat) for plotting
    overlaydata(:,:) = overlaydata_raw(:,:);
    overlaydata(:,2) = sin(pi*overlaydata_raw(:,2)/180.0);
    % grid (and average per cell) data if requested
    % NOTE: data vector length is re-calculated and the value of nmax reset
    if (data_ijk_mean == 'y')
        overlaydata_ij_old(:,:) = overlaydata_ij(:,:);
        overlaydata_ij(:,:) = [];
        overlaydata(:,:)    = [];
        m=0;
        for i = 1:imax,
            for j = 1:jmax,
                if (~isnan(zm(j,i)))
                    samecell_locations = find((overlaydata_ij_old(:,1)==i)&(overlaydata_ij_old(:,2)==j));
                    samecell_n = size(samecell_locations);
                    if (samecell_n(1) > 0)
                        m=m+1;
                        samecell_mean = mean(overlaydata_ij_old(samecell_locations,3));
                        overlaydata_ij(m,:) = [i j samecell_mean];
                        overlaydata(m,1) = grid_lon_origin + 360.0*(overlaydata_ij(m,1) - 0.5)/jmax;
                        overlaydata(m,2) = 2.0*(overlaydata_ij(m,2) - 0.5)/jmax - 1.0;
                        overlaydata(m,3) = samecell_mean;
                    end
                end
            end
        end
        nmax=m;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** TAYLOR DIAGRAM **************************************************** %
% *********************************************************************** %
% calculate stats needed for Taylor Diagram (and plot it!)
% *********************************************************************** %
%
% *** DISCRETE DATA ***************************************************** %
%
% NOTE: no scale transformatoin has been appplied
%       to either gridded or % overlay data
% NOTE: valid only for data on a single depth level
if (~isempty(overlaydataid) && ((data_only == 'n') || (data_anomoly == 'y')))
    % set overlay data vector
    data_vector_1 = overlaydata(:,3);
    % populate the gridded dataset vector with values corresponding to
    % the overlay data locations
    % NOTE: !!! data is (j,i) !!! (=> swap i and j)
    % NOTE: re-orientate data_vector_2 to match data_vector_1
    for n = 1:nmax,
        data_vector_2(n) = data(overlaydata_ij(n,2),overlaydata_ij(n,1));
    end
    data_vector_2 = data_vector_2';
    % filter data
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    if (data_stats == 'y')
        % calculate stats
        STATM = calc_allstats(data_vector_1,data_vector_2);
        % plot Taylor diagram
        taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
        print('-depsc2', [filename, '_TaylorDiagram.', str_date, '.eps']);
        %%%% plot Target diagram
        %%%targetdiag_vargout = plot_target(STATM(7,1:2),STATM(8,1:2),'r',1.0,[],[]);
        %%%print('-depsc2', [filename, '_TargetDiagram.', str_date, '.eps']);
    end
else
    STATM = [];
end
%
% *** SAVE STATS DATA *************************************************** %
%
% STATM(1,:) = MEAN;
% STATM(2,:) = STD;
% STATM(3,:) = RMSD;
% STATM(4,:) = CORRELATIONS;
% STATM(5,:) = N;
% STATM(6,:) = TOTAL RMSD;
% STATM(7,:) = STANDARD DEVIATION NORMALISED RMSD;
% STATM(8,:) = NORMALISED BIAS;
% STATM(9,:) = R2;
if (data_stats == 'y')
    if (~isempty(overlaydataid) && ((data_only == 'n') || (data_anomoly == 'y')))
        fid = fopen([filename '_STATS' '.dat'], 'wt');
        fprintf(fid, '\n');
        fprintf(fid, '=== STATS SUMMARY ===');
        fprintf(fid, '\n');
        fprintf(fid, 'Number of data points, N                           : %8.6e \n', STATM(5,2));
        fprintf(fid, '\n');
        fprintf(fid, 'Mean                                               : %8.6e \n', STATM(1,2));
        fprintf(fid, 'R2                                                 : %8.6e \n', STATM(9,2));
        fprintf(fid, '\n');
        fprintf(fid, 'Standard Deviation (scaled by N)                   : %8.6e \n', STATM(2,2));
        fprintf(fid, 'Root Mean Square Difference (scaled by N)          : %8.6e \n', STATM(3,2));
        fprintf(fid, 'TOTAL Root Mean Square Difference                  : %8.6e \n', STATM(6,2));
        fprintf(fid, 'Correlation                                        : %8.6e \n', STATM(4,2));
        fprintf(fid, 'STANDARD DEVIATION NORMALISED RMSD                 : %8.6e \n', STATM(7,2));
        fprintf(fid, 'NORMALISED BIAS                                    : %8.6e \n', STATM(8,2));
        fclose(fid);
    end
end
%
% *** SAVE EQUIVALENT MODEL DATA **************************************** %
%
% save model data at the data locations
if (~isempty(overlaydataid) && (data_only == 'n'))
    fid = fopen([filename '_MODELPOINTS', '.', str_date '.dat'], 'wt');
    fprintf(fid, '%% Model value at data locations');
    fprintf(fid, '\n');
    fprintf(fid, '%% Format: i, j, lon, lat, model value, data value, label');
    fprintf(fid, '\n');
    for n = 1:nmax,
        fprintf(fid, '%d %d %8.3f %8.3f %8.6e %8.6e %s \n', int16(overlaydata_ij(n,1)), int16(overlaydata_ij(n,2)), overlaydata(n,1), 180.0*asin(overlaydata(n,2))/pi, data_vector_2(n), overlaydata(n,3), overlaylabel(n,:));
    end
    fclose(fid);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** ANOMOLY PLOTTING DATA ADJUSTMENTS ********************************* %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % calculate molde-data anomoly
    if (data_anomoly == 'y')
        overlaydata(:,3) = data_vector_2(:) - data_vector_1(:);
    end
    % redefine model grid location values so as to all appear white
    if (data_only == 'y')
        zm = zeros(size(zm(:,:)));
        zm(find(zm(:,:) == 0)) = NaN;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** TRANSFORM LON GRID ************************************************ %
% *********************************************************************** %
%
% extend gridded data in +/- longitude
lon_start = min(min(lonw));
i_start = round((lon_min-(lon_start-360.0))/(360.0/jmax)) + 1;
xm_ex = [xm - 360.0 xm + 000.0 xm + 360.0];
ym_ex = [ym + 000.0 ym + 000.0 ym + 000.0];
zm_ex = [zm zm zm];
topo_ex = [topo topo topo];
lonm_ex = [lonm - 360.0 lonm + 000.0 lonm + 360.0];
lone_ex = [lone - 360.0 lone + 000.0 lone + 360.0];
lonw_ex = [lonw - 360.0 lonw + 000.0 lonw + 360.0];
layb_ex = [layb layb layb];
% shorten to conform to desired lon start value
xm = xm_ex(:,i_start:i_start+imax-1);
ym = ym_ex(:,i_start:i_start+imax-1);
zm = zm_ex(:,i_start:i_start+imax-1);
topo = topo_ex(:,i_start:i_start+imax-1);
lonm = lonm_ex(:,i_start:i_start+imax-1);
lone = lone_ex(:,i_start:i_start+imax-1);
lonw = lonw_ex(:,i_start:i_start+imax-1);
layb = layb_ex(:,i_start:i_start+imax-1);
if ~isempty(overlaydataid)
    % force discrete data to lie within longitude plotting axis
    % (lon_min to lon_min + 360)
    for n = 1:nmax,
        if (overlaydata(n,1) < lon_min)
            overlaydata(n,1) = overlaydata(n,1) + 360.0;
        end
        if (overlaydata(n,1) > (lon_min + 360.0))
            overlaydata(n,1) = overlaydata(n,1) - 360.0;
        end
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT MAIN FIGURE ************************************************** %
% *********************************************************************** %
%
% *** CONFIGURE AND CREATE PLOTTING WINDOW ****************************** %
%
% create figure
scrsz = get(0,'ScreenSize');
figure('Position',[((1 - plot_dscrsz)/2)*plot_dscrsz*scrsz(3) (1 - plot_dscrsz)*plot_dscrsz*scrsz(4) plot_dscrsz*scrsz(3) 0.60*plot_dscrsz*scrsz(4)])
clf;
% define plotting regions
if (plot_format_old == 'y')
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.10 0.05 0.65 0.90]);
    fh(3) = axes('Position',[0.80 0.27 0.20 0.46],'Visible','off');
else
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.15 0.15 0.65 0.70]);
    fh(3) = axes('Position',[0.75 0.15 0.15 0.70],'Visible','off');
end
% define colormap
cmap = make_cmap(colorbar_name,con_n+2);
if (colorbar_inv == 'y'), cmap = flipdim(cmap,1); end,
colormap(cmap);
% date-stamp plot Dom: could give units here:
set(gcf,'CurrentAxes',fh(1));
if (plot_format_old == 'y')
    text(0.95,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
else
    text(0.85,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
end
%
% *** SET PLOT SCALE **************************************************** %
%
% set minimum contour value
if exist('con_min','var') == 0
    con_min = min(min(zm));
end
% set maximum contour value
if exist('con_max','var') == 0
    con_max = max(max(zm));
end
% ensure min and max are not identical ...
if con_min == con_max
    if con_max == 0.0
        con_max = 1.0;
    else
        con_min = (1.00/1.01)*con_min;
        con_max = (1.01)*con_max;
    end
end
% if min > max, then reverse min and max
if con_min > con_max
    con_min_TEMP = con_min;
    con_max_TEMP = con_max;
    con_min = con_max_TEMP;
    con_max = con_min_TEMP;
end
%
% *** CREATE MAIN PLOT ************************************************** %
%
set(gcf,'CurrentAxes',fh(2));
hold on;
% set color and lat/lon axes and labels
caxis([con_min-(con_max-con_min)/con_n con_max]);
set(gca,'PlotBoxAspectRatio',[1.0 plot_xy_scaling*0.5 1.0]);
axis([lon_min lon_max lat_min lat_max]);
if plot_global,
    axis([lon_min lon_max lat_min lat_max]);
    set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',[lon_min:plot_lon_delta:lon_max]);
    if (plot_equallat == 'n'),
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
    else
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-90.0 -60.0 -30.0 0 30.0 60.0 90.0], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
    end
else
    axis([plot_lon_min plot_lon_max plot_lat_min plot_lat_max]);
    set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',[plot_lon_min plot_lon_max]);
    if (plot_equallat == 'n'),
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[plot_lat_min plot_lat_max], 'YTickLabel',{num2str(180*asin(plot_lat_min)/pi);num2str(180*asin(plot_lat_max)/pi)});
    else
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[plot_lat_min plot_lat_max], 'YTickLabel',{num2str(plot_lat_min);num2str(plot_lat_max)});
    end
end
set(gca,'TickDir','out');
if isempty(plot_title)
    plot_title = [' ',strrep(dataid_1,'_',' ')];
    plot_title(find(plot_title(:)=='_')) = '-';
end
title(plot_title,'FontSize',12);
% draw filled rectangles: Dom this is were results are drawn!
for i = 1:imax,
    for j = 1:jmax,
        if topo(j,i) > layb(j,i)
            h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],color_g);
            set(h,'EdgeColor',color_g);
        else
            if (isnan(zm(j,i)))
                if exist('data_nancolor')
                    h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],data_nancolor);
                    set(h,'EdgeColor',data_nancolor);
                else
                    h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],[1 1 1]);
                    set(h,'EdgeColor',[1 1 1]);
                end
            else    % Dom: create color depending on value
                col = 1 + round(0.5+con_n*(zm(j,i)-con_min)/(con_max-con_min));
                if col < 1, col = 1; end
                if col > con_n+2, col = con_n+2; end
                h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],cmap(col,:));
                set(h,'EdgeColor',cmap(col,:));
            end
        end
    end
end
%
% *** PLOT CONTINENTAL OUTLINE ****************************************** %
%
% draw continental outline
for j = 1:jmax,
    for i = 1:imax-1,
        if topo(j,i) > layb(j,i)
            if topo(j,i+1) <= layb(j,i+1)
                h = plot([lone(j,i) lone(j,i)],[lats(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
    for i = 2:imax,
        if topo(j,i) > layb(j,i)
            if topo(j,i-1) <= layb(j,i-1)
                h = plot([lonw(j,i) lonw(j,i)],[lats(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
for i = 1:imax,
    for j = 1:jmax-1,
        if topo(j,i) > layb(j,i)
            if topo(j+1,i) <= layb(j+1,i)
                h = plot([lonw(j,i) lone(j,i)],[latn(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
    for j = 2:jmax,
        if topo(j,i) > layb(j,i)
            if topo(j-1,i) <= layb(j-1,i)
                h = plot([lonw(j,i) lone(j,i)],[lats(j,i) lats(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
%
% *** OVERLAY CONTOURS ************************************************** %
%
% plot contours
if (contour_plot == 'y') && (data_only == 'n'),
    v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
    [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
    set(h,'LineWidth',0.25);
    v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
    [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k');
    set(h,'LineWidth',0.5);
    if data_log10 == 'y'
        %%%%%%%%
    elseif contour_label == 'y'
        clabel(C,h);
    end
end
%
% *** OVERLAY DATA ****************************************************** %
%
if ~isempty(overlaydataid)
    % set uniform marker shape and color
    if (data_shapecol == 'n'),
        for n = 1:nmax,
            overlaydata_shape(n) = 'o';
            overlaydata_fcol(n) = data_sitecolor;
            if exist('data_linecolor')
                overlaydata_ecol(n) = data_linecolor;
            else
                overlaydata_ecol(n) = data_sitecolor;
            end
        end
    end
    % plot overlay data
    if (data_siteonly == 'n')
        scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata(:,3)/data_scale,overlaydata_shape(n),'Filled','LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
    else
        if (overlaydata_fcol(n) == '-'),
            scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
        else
            scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n),'MarkerFaceColor',overlaydata_fcol(n));
        end
    end
    if (data_sitelabel == 'y'),
        text(overlaydata(:,1)+(data_size/30),overlaydata(:,2)+(data_size/1200),overlaylabel(:,:),'FontSize',data_fontsz,'Color',data_sitecolor);    
    end
end
%
% *** PLOT BORDER ******************************************************* %
%
% draw plot border
h = plot([lon_min lon_max],[lat_min lat_min],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_min lon_max],[lat_max lat_max],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_min lon_min],[lat_min lat_max],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_max lon_max],[lat_min lat_max],'k-');
set(h,'LineWidth',1.0);
%
hold off;
%
% *** CREATE COLOR BAR ************************************************** %
%
if (~((data_only == 'y') && (data_siteonly == 'y')))
    %
    set(gcf,'CurrentAxes',fh(3));
    hold on;
    %
    set(gca,'XTick',[],'YTick',[]);
    axis([0 1 0 con_n+2]);
    % draw and label color bar rectangles
    % draw and label start triangle
    c = 1;
    h = fill([0.1 0.2 0.3],[c c-1.0 c],cmap(c,:));
    if isempty(contour_file),
        str = [num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
    else
        str = num2str(contour_data(c));
    end
    textsize = 2+round(80/con_n);
    if textsize > 10, textsize = 10; end
    text(0.40,c,str,'FontName','Arial','FontSize',textsize);
    set(h,'LineWidth',0.5);
    set(h,'EdgeColor','k');
    % draw and label bars
    for c = 2:con_n+1,
        h = fill([0.1 0.1 0.3 0.3],[c-1.0 c c c-1.0],cmap(c,:));
        if isempty(contour_file),
            str = [num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
        else
            str = num2str(contour_data(c));
        end
        textsize = 2+round(80/con_n);
        if textsize > 10, textsize = 10; end
        text(0.40,c,str,'FontName','Arial','FontSize',textsize);
        set(h,'LineWidth',0.5);
        set(h,'EdgeColor','k');
    end
    % draw end triangle
    c = con_n+2;
    h = fill([0.1 0.2 0.3],[c-1.0 c c-1.0],cmap(c,:));
    set(h,'LineWidth',0.5);
    set(h,'EdgeColor','k');
    %
    hold off;
    %
    hold off;
    %
end
%
% *** PRINT PLOT ******************************************************** %
%
set(gcf,'CurrentAxes',fh(1));
if (plot_format_old == 'y')
    print('-dpsc2', [filename '.' str_date '_20_colors.ps']);
else
    switch plot_format
        case 'png'
            export_fig([filename '.' str_date '.png'], '-png', '-r150', '-nocrop');
        case 'pngT'
            export_fig([filename '.' str_date '.png'], '-png', '-r150', '-nocrop', '-transparent');
        case 'jpg'
            export_fig([filename '.' str_date '.jpg'], '-jpg', '-r150', '-nocrop');
        otherwise
            export_fig([filename '.' str_date '.eps'], '-eps', '-nocrop');
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SECONDARY FIGURES ************************************************* %
% *********************************************************************** %
%
% % % if (plot_secondary == 'y'),
% % %     %
% % %     % *** PLOT FIGURE (cross-plot) ************************************** %
% % %     %
% % %     if ( ~isempty(dataid_2) || (~isempty(overlaydataid) && (data_only == 'n')) ),
% % %         %
% % %         if ~isempty(dataid_2),
% % %             loc_x_data = data_vector_1;
% % %             loc_y_data = data_vector_2;
% % %             loc_x_label = [strrep(dataid_1,'_','-')];
% % %             loc_y_label = [strrep(dataid_2,'_','-')];
% % %         elseif ~isempty(overlaydataid)
% % %             loc_x_data = data_vector_1;
% % %             loc_y_data = data_vector_2;
% % %             loc_x_label = [strrep(overlaydataid,'_','-')];
% % %             loc_y_label = [strrep(dataid_1,'_','-')];
% % %         end
% % %         %
% % %         plot_crossplotc(loc_x_data,loc_y_data,[],loc_x_label,loc_y_label,'',POPT,[filename '.CROSSPLOT']);
% % %         %
% % %     end
% % %     %
% % %     % *** SAVE DATA (cross-plot relationships) ************************** %
% % %     %
% % %     if ( ~isempty(dataid_2) || (~isempty(overlaydataid) && (data_only == 'n')) ), fprint_1Dn_d([loc_x_data loc_y_data],[filename '.CROSSPLOT.', str_date, '.res']); end
% % %     %
% % % end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% close netCDF files
netcdf.close(ncid_1);
if ~isempty(exp_2)
    netcdf.close(ncid_2);
end
%
% *********************************************************************** %
